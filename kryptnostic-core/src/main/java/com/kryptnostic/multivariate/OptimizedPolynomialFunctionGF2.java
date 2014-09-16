package com.kryptnostic.multivariate;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.Executors;

import org.apache.commons.lang3.tuple.Pair;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import cern.colt.bitvector.BitVector;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.google.common.util.concurrent.ListeningExecutorService;
import com.google.common.util.concurrent.MoreExecutors;
import com.kryptnostic.linear.BitUtils;
import com.kryptnostic.linear.EnhancedBitMatrix;
import com.kryptnostic.multivariate.gf2.Monomial;
import com.kryptnostic.multivariate.gf2.SimplePolynomialFunction;

public class OptimizedPolynomialFunctionGF2 extends BasePolynomialFunction {
    private static final Logger logger = LoggerFactory.getLogger(OptimizedPolynomialFunctionGF2.class);
    protected static final int CONCURRENCY_LEVEL = Runtime.getRuntime().availableProcessors() - 1;
    protected static final ListeningExecutorService executor = MoreExecutors.listeningDecorator(Executors
            .newFixedThreadPool(CONCURRENCY_LEVEL));

    public OptimizedPolynomialFunctionGF2(int inputLength, int outputLength, Monomial[] monomials,
            BitVector[] contributions) {
        super(inputLength, outputLength, monomials, contributions);
    }

    public BitVector apply(final BitVector input) {

        final CountDownLatch latch = new CountDownLatch(CONCURRENCY_LEVEL);

        final BitVector result = new BitVector(outputLength);
        int blocks = ( monomials.length / CONCURRENCY_LEVEL );
        int leftover = monomials.length % CONCURRENCY_LEVEL;

        for (int i = 0; i < CONCURRENCY_LEVEL; i++) {
            final int fromIndex = i * blocks;
            int targetIndex = fromIndex + blocks;
            if (leftover != 0 && i == CONCURRENCY_LEVEL - 1) {
                targetIndex += leftover;
            }
            final int toIndex = targetIndex;

            Runnable r = new Runnable() {
                @Override
                public void run() {
                    BitVector intermediary = new BitVector(outputLength);
                    for (int i = fromIndex; i < toIndex; ++i) {
                        Monomial term = monomials[i];
                        if (term.eval(input)) {
                            intermediary.xor(contributions[i]);
                        }
                    }
                    synchronized (result) {
                        result.xor(intermediary);
                    }
                    latch.countDown();
                }
            };
            executor.execute(r);
        }
        try {
            latch.await();
        } catch (InterruptedException e) {
            logger.error("Concurrent apply() latch interrupted.");
        }
        return result;
    }

    /**
     * Composes an outer function with the inner function.
     * 
     */
    @Override
    public SimplePolynomialFunction compose(SimplePolynomialFunction inner) {
        Preconditions.checkArgument(inputLength == inner.getOutputLength(),
                "Input length of outer function must match output length of inner function it is being composed with");

        if (this.getMaximumMonomialOrder() == 2 && inner.getMaximumMonomialOrder() == 1) {
            return composeQuadratic(inner);
        }
        ComposePreProcessResults prereqs = preProcessCompose(inner);

        BitVector[] results = expandOuterMonomials(prereqs.monomialsList, prereqs.innerRows, prereqs.indices);

        return postProcessCompose(prereqs.monomialsList, prereqs.indices, results, inner);
    }

    /**
     * Factors outer terms by first set bit, so that product is linear.
     * 
     * @return SimplePolynomialFunction
     */
    private SimplePolynomialFunction composeQuadratic(SimplePolynomialFunction inner) {
        // group outer monomial-contribution pairs by common factor
        Pair<List<Integer>, List<List<Integer>>> subfunctions = factorSubFunctions();
        Pair<List<Monomial>, ConcurrentMap<Monomial, Integer>> mappedMonomials = getMonomialMap(inner);

        // expand each of these subfunctions
        List<Integer> factors = subfunctions.getLeft();
        List<List<Integer>> outerSubFunctions = subfunctions.getRight();
        List<BitVector[]> newContributions = Lists.newArrayList();
        for (int i = 0; i < factors.size(); i++) {
            // TODO handle constant outer monomial separately -- it cannot be a factor
            BitVector[] newContribution = expandSubfunction(factors.get(i), outerSubFunctions.get(i), inner,
                    mappedMonomials.getLeft(), mappedMonomials.getRight());
            newContributions.add(newContribution);
        }

        return composeReduce(newContributions, mappedMonomials.getLeft(), inner);
    }

    private Pair<List<Monomial>, ConcurrentMap<Monomial, Integer>> getMonomialMap(SimplePolynomialFunction inner) {
        List<Monomial> mList = Lists.newArrayList(inner.getMonomials());
        ConcurrentMap<Monomial, Integer> indices = Maps.newConcurrentMap();
        for (int i = 0; i < mList.size(); ++i) {
            indices.put(mList.get(i), i);
        }
        return Pair.of(mList, indices);
    }

    /**
     * @return function containing collected monomials and contributions from each subfunction expansion
     */
    private SimplePolynomialFunction composeReduce(List<BitVector[]> bucketedContributions, List<Monomial> mList,
            SimplePolynomialFunction inner) {
        // TODO rewrite lengths of contributions, filter contributions and monomials (?)
        Monomial[] monomials = mList.toArray(new Monomial[0]);
        SimplePolynomialFunction composed = new BasePolynomialFunction(inner.getInputLength(), outputLength, monomials,
                contributions);
        return composed;
    }

    /**
     * @return function representing expanded outer monomials grouped by common factor.
     */
    private BitVector[] expandSubfunction(Integer commonFactor, List<Integer> monomialIndices,
            SimplePolynomialFunction inner, List<Monomial> mList, ConcurrentMap<Monomial, Integer> indices) {
        // Select outer monomials and contributions by indices
        EnhancedBitMatrix outerContributions = getOuterContributions(monomialIndices);
        EnhancedBitMatrix innerContributions = getInnerContributions(commonFactor, monomialIndices,
                inner.getContributions());
        // inner x outer contributions
        EnhancedBitMatrix linearContributions = outerContributions.multiply(innerContributions);
        // compute product of common-factor with expanded linear contributions
        BitVector commonFactorContribution = inner.getContributions()[commonFactor];
        List<BitVector> linearContributionRows = linearContributions.getRows();
        BitVector[] newContributions = new BitVector[linearContributionRows.size()];
        for (int i = 0; i < linearContributionRows.size(); i++) {
            newContributions[i] = product(commonFactorContribution, linearContributionRows.get(i), mList, indices);
        }
        return newContributions;
    }

    /**
     * @param bitVectors
     * @return EnhancedBitMatrix of inner contributions corresponding to the non-shared factors
     * 
     *         m' x Iol, where m' is the number of distinct, non-shared factors.
     */
    private EnhancedBitMatrix getInnerContributions(Integer commonFactor, List<Integer> monomialIndices,
            BitVector[] innerContributions) {
        Set<Integer> factors = Sets.newHashSet();
        for (Integer index : monomialIndices) {
            Monomial m = monomials[index];
            List<Integer> secondaryFactors = BitUtils.assertedIndices(m, commonFactor + 1);
            factors.addAll(secondaryFactors);
        }
        // Get inner contributions in row form -- outputlength x numterms
        EnhancedBitMatrix contributionRows = new EnhancedBitMatrix(Lists.newArrayList(innerContributions));
        contributionRows = contributionRows.tranpose();
        List<BitVector> innerRows = contributionRows.getRows();
        // Select the rows corresponding to the outer monomial factor for inclusion in result matrix
        List<BitVector> selectedRows = Lists.newArrayList();
        for (Integer factor : factors) {
            selectedRows.add(innerRows.get(factor));
        }
        EnhancedBitMatrix result = new EnhancedBitMatrix(selectedRows);
        return result.tranpose();
    }

    /**
     * @return EnhancedBitMatrix of the outer contributions corresponding to the subset of outer monomials which share a
     *         common factor.
     * 
     *         Ool x m, where m is the number of monomials sharing a factor
     */
    private EnhancedBitMatrix getOuterContributions(List<Integer> monomialIndices) {
        List<BitVector> contributionRows = Lists.newArrayList();
        for (Integer index : monomialIndices) {
            contributionRows.add(contributions[index].copy());
        }
        return new EnhancedBitMatrix(contributionRows);
    }

    /**
     * @return Pair of an array of outer monomial factors and an array of Sets of indices of the outer monomials sharing
     *         these factors.
     */
    private Pair<List<Integer>, List<List<Integer>>> factorSubFunctions() {
        List<Integer> factors = Lists.newArrayList();
        List<List<Integer>> bucketIndices = Lists.newArrayList();
        Map<Integer, Integer> factorIndices = Maps.newHashMap();

        for (int i = 0; i < monomials.length; i++) {
            Monomial m = monomials[i];
            if (m.isZero()) {
                continue;
                // TODO handle constant monomials
            }
            int factor = BitUtils.first(m);

            List<Integer> bucket;
            Integer index = factorIndices.get(factor);
            if (index == null) {
                index = factors.size();
                factors.add(index);
                factorIndices.put(factor, index);
                bucket = Lists.newArrayList();
                bucketIndices.add(bucket);
            } else {
                bucket = bucketIndices.get(index);
            }
            bucket.add(i);

        }
        return Pair.of(factors, bucketIndices);
    }

    @Override
    protected BitVector[] expandOuterMonomials(final List<Monomial> mList, final BitVector[] innerRows,
            final ConcurrentMap<Monomial, Integer> indices) {
        final CountDownLatch latch = new CountDownLatch(CONCURRENCY_LEVEL);
        final BitVector[] results = new BitVector[monomials.length];
        int blocks = monomials.length / CONCURRENCY_LEVEL;
        int leftover = monomials.length % CONCURRENCY_LEVEL;

        for (int i = 0; i < CONCURRENCY_LEVEL; i++) {
            final int fromIndex = i * blocks;
            int targetIndex = fromIndex + blocks;
            if (leftover != 0 && i == CONCURRENCY_LEVEL - 1) {
                targetIndex += leftover;
            }
            final int toIndex = targetIndex;

            executor.execute(new Runnable() {

                @Override
                public void run() {
                    for (int j = fromIndex; j < toIndex; j++) {
                        Monomial outerMonomial = monomials[j];
                        BitVector newContributions = null;
                        if (outerMonomial.isZero()) {
                            newContributions = new BitVector(mList.size());
                        } else {
                            for (int i = Long.numberOfTrailingZeros(outerMonomial.elements()[0]); i < outerMonomial
                                    .size(); ++i) {
                                if (outerMonomial.get(i)) {
                                    if (newContributions == null) {
                                        newContributions = innerRows[i];
                                    } else {
                                        newContributions = product(newContributions, innerRows[i], mList, indices);
                                    }
                                }
                            }
                        }
                        results[j] = newContributions;
                    }
                    latch.countDown();
                }
            });

        }

        try {
            latch.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        return results;
    }
}

package com.kryptnostic.multivariate;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.Executors;

import org.apache.commons.lang3.tuple.Pair;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import cern.colt.bitvector.BitVector;

import com.google.common.base.Optional;
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
import com.kryptnostic.multivariate.parameterization.ParameterizedPolynomialFunctionGF2;

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
        Pair<SortedMap<Integer, Set<Integer>>, BitVector[]> factoringResults = factorOuterMonomials();
        SortedMap<Integer, Set<Integer>> factors = factoringResults.getLeft();
        BitVector[] reducedOuterContributions = factoringResults.getRight();
        // get prereqs
        ComposePreProcessResults prereqs = preProcessCompose(inner);
        // expand inner factors
        BitVector[] reducedInnerContributions = reduceInnerContributions(factors, prereqs.innerRows);
        BitVector[] results = expandFactoredOuterMonomials(factors.keySet(), reducedInnerContributions, prereqs.monomialsList,
                prereqs.innerRows, prereqs.indices);

        return postProcessFactoredCompose(prereqs.monomialsList, prereqs.indices, results, inner, reducedOuterContributions);
    }

    /**
     * Constructs an array of the inner contributions
     * 
     * @return BitVector[]
     */
    private BitVector[] reduceInnerContributions(SortedMap<Integer, Set<Integer>> factors, BitVector[] innerRows) {
        BitVector[] reducedContributions = new BitVector[factors.keySet().size()];
        int i = 0;
        for (Integer firstFactor : factors.keySet()) {
            Set<Integer> secondFactors = factors.get(firstFactor);
            // Sum contributions with shared factor.
            BitVector innerContributionSum = new BitVector(innerRows[0].size());
            for (Integer secondFactor : secondFactors) {
                if (secondFactor == -1) {
                    innerContributionSum.not(); // second factor is constant, so invert the sum
                } else {
                    innerContributionSum.xor(innerRows[secondFactor]);
                }
            }
            reducedContributions[i] = innerContributionSum;
            i++;
        }
        return reducedContributions;
    }

    /**
     * Buckets outer monomials by first set bit, and returns a Map of first set bit to a set of second set bits in the
     * monomials array.
     * 
     * @return Map<Integer,Set<Integer>>
     */
    private Pair<SortedMap<Integer, Set<Integer>>, BitVector[]> factorOuterMonomials() {
        SortedMap<Integer, Set<Integer>> factors = Maps.newTreeMap();
        Map<Integer, BitVector> newContributions = Maps.newHashMap();
        for (int i = 0; i < monomials.length; i++) {
            Monomial m = monomials[i].clone();
            if (m.isZero()) {
                continue;
            }
            
            Integer firstSetBit = BitUtils.first(m);
            Integer secondSetBit = -1;
            m.clear(firstSetBit);
            if (!m.isZero()) {
                secondSetBit = BitUtils.first(m);
            }

            Set<Integer> secondSetBits = factors.get(firstSetBit);
            if (secondSetBits == null) {
                secondSetBits = Sets.newHashSet();
                
            } else {
                
            }
            
            BitVector reducedOuterContribution = newContributions.get(firstSetBit);
            if (reducedOuterContribution == null) {
                reducedOuterContribution = contributions[i];
                newContributions.put(firstSetBit, reducedOuterContribution);
            } else {
                reducedOuterContribution.xor(contributions[i]);
                newContributions.put(firstSetBit, reducedOuterContribution);
            }
            
            secondSetBits.add(secondSetBit);
            factors.put(firstSetBit, secondSetBits);
        }
        
        return Pair.of(factors, newContributions.values().toArray(new BitVector[newContributions.values().size()]));
    }

    // TODO make concurrent
    protected BitVector[] expandFactoredOuterMonomials(Iterable<Integer> outerFactors, BitVector[] reducedRows,
            final List<Monomial> mList, final BitVector[] innerRows, final ConcurrentMap<Monomial, Integer> indices) {
        List<BitVector> results = Lists.newArrayList();
        int i = 0;
        for (Integer factor : outerFactors) {
            BitVector newContribution = product(innerRows[factor], reducedRows[i], mList, indices);
            results.add(newContribution);
            i++;
        }
        return results.toArray(new BitVector[results.size()]);
    }
    
    private SimplePolynomialFunction postProcessFactoredCompose(List<Monomial> mList, Map<Monomial, Integer> indices,
            BitVector[] results, SimplePolynomialFunction inner, BitVector[] reducedOuterContributions) {
                
        Optional<Integer> constantOuterMonomialIndex = Optional.absent();
        Optional<Integer> constantInnerMonomialIndex = Optional.fromNullable(indices.get(Monomial
                .constantMonomial(inner.getInputLength())));
        
        // Now lets fix the contributions so they're all the same length.
        for (int i = 0; i < results.length; ++i) {
            BitVector contribution = results[i];
            if (contribution.size() != mList.size()) {
                contribution.setSize(mList.size());
            }
        }
        /*
         * Calculate resulting set of contributions in terms of the new monomial basis.
         * So for every outer monomial factor, xor all of the contributions of the outer monomials that share that factor,
         * then calculate new outer contributions using these condensed contributions and the results.
         */
        BitVector[] outputContributions = new BitVector[outputLength];

        for (int row = 0; row < outputLength; row++) {
            outputContributions[row] = new BitVector(mList.size());
            for (int i = 0; i < reducedOuterContributions.length; i++) {
                if (reducedOuterContributions[i].get(row)) {
                    if (monomials[i].isZero()) {
                        constantOuterMonomialIndex = Optional.of(i);
                    } else {
                        outputContributions[row].xor(results[i]);
                    }
                }
            }
        }

        /*
         * After we have computed the contributions in terms of the new monomial basis we transform from row to column
         * form of contributions to match up with each monomial in mList
         */
        List<BitVector> unfilteredContributions = Lists.newArrayList(outputContributions);
        EnhancedBitMatrix.transpose(unfilteredContributions, mList.size());

        /*
         * If the outer monomial has constant terms and the unfiltered contributions have a constant term, than we xor
         * them together to get the overall constant contributions.
         */

        if (constantOuterMonomialIndex.isPresent()) {
            if (constantInnerMonomialIndex.isPresent()) {
                unfilteredContributions.get(constantInnerMonomialIndex.get()).xor(
                        contributions[constantOuterMonomialIndex.get()]);
            } else {
                // Don't use the outer monomial directly since it maybe the wrong size.
                // mList.add( monomials[ constantOuterMonomialIndex.get() ] );
                mList.add(Monomial.constantMonomial(inner.getInputLength()));
                unfilteredContributions.add(contributions[constantOuterMonomialIndex.get()]);
            }
        }

        /*
         * Now we filter out any monomials, which have nil contributions.
         */

        List<BitVector> filteredContributions = Lists.newArrayListWithCapacity(unfilteredContributions.size());
        List<BitVector> filteredMonomials = Lists.newArrayListWithCapacity(mList.size());
        for (int i = 0; i < mList.size(); ++i) {
            BitVector contrib = unfilteredContributions.get(i);
            if (notNilContributionPredicate.apply(contrib)) {
                filteredContributions.add(contrib);
                filteredMonomials.add(mList.get(i));
            }
        }

        if (inner.isParameterized()) {
            ParameterizedPolynomialFunctionGF2 ppf = (ParameterizedPolynomialFunctionGF2) inner;
            return new ParameterizedPolynomialFunctionGF2(inner.getInputLength(), outputLength,
                    filteredMonomials.toArray(new Monomial[0]), filteredContributions.toArray(new BitVector[0]),
                    ppf.getPipelines());
        }

        return new BasePolynomialFunction(inner.getInputLength(), outputLength,
                filteredMonomials.toArray(new Monomial[0]), filteredContributions.toArray(new BitVector[0]));
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

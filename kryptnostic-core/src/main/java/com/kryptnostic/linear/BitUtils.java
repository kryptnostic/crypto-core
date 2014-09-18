package com.kryptnostic.linear;

import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.Set;

import cern.colt.bitvector.BitVector;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.kryptnostic.multivariate.MultivariateUtils;

public final class BitUtils {
    private BitUtils() {
    }

    // TODO:Re-enable seeding.
    private static final Random r = new Random(0);// System.currentTimeMillis() );

    public static long parity(long l) {
        l ^= l >> 32;
        l ^= l >> 16;
        l ^= l >> 8;
        l ^= l >> 4;
        l ^= l >> 2;
        l ^= l >> 1;
        return l & 1L;
    }

    public static int getFirstSetBit(BitVector v) {
        // TODO: Optimize
        for (int i = 0; i < v.size(); ++i) {
            if (v.get(i)) {
                return i;
            }
        }
        return -1;
    }

    public static BitVector randomVector(int length) {
        return MultivariateUtils.randomVector(length);
    }

    /**
     * @param from
     *            , the index in the backing long array to start from
     * @param to
     *            , the last index in the backing long array to copy
     * @return BitVector
     */
    // TODO consider refactoring this so as not to reach into implementation.
    public static BitVector subVector(BitVector v, int from, int to) {
        return new BitVector(Arrays.copyOfRange(v.elements(), from, to), ( to - from ) << 6);
    }
    
    /**
     * Constructs a new BitVector from the bit range specified, where
     * indices refer to the indices the given BitVector.
     * @return BitVector
     */
    public static BitVector newFromTo(BitVector v, int fromIndex, int toIndex) {
        Preconditions.checkArgument(toIndex > fromIndex, "toIndex must be larger than fromIndex");
        Preconditions.checkArgument(toIndex < v.size(), "toIndex must be less than size of original vector.");
        BitVector result = new BitVector(toIndex + 1 - fromIndex);
        for (int i = fromIndex; i <= toIndex; i++) {
            if (v.get(i)) {
                result.set(i - fromIndex);
            }
        }
        return result;
    }

    public static BitVector randomVector(int length, int desiredHammingWeight) {
        BitVector v = new BitVector(length);
        /*
         * In theory optimized popcnt instruction is going to be faster than bit twiddling to check individual bits.
         */
        while (v.cardinality() < desiredHammingWeight) {
            v.set(r.nextInt(length));
        }
        return v;
    }

    public static BitVector extend(BitVector v, int newSize) {
        Preconditions.checkArgument(v.size() <= newSize, "New size must be greater than input vector size.");
        BitVector result = v.copy();
        result.setSize(newSize);
        return result;
    }

    /**
     * Given a mapping from old indices to new indices, creates a new BitVector ordered by this mapping.
     * 
     * @return BitVector
     */
    public static BitVector extendAndOrder(BitVector v, int[] mapping, int newSize) {
        Preconditions.checkArgument(mapping.length == v.size(), "Must map exactly every bit in the vector to a new index.");
        Set<Integer> elements = Sets.newHashSet();
        BitVector sorted = new BitVector(newSize);
        for (int i = 0; i < v.size(); i++) {
            int index = mapping[i];
            Preconditions.checkArgument(!elements.contains(index), "Cannot map two variables to the same index.");
            elements.add(index);
            Preconditions.checkArgument(index < newSize, "Cannot map to index greater than new size.");
            if (v.get(i)) {
                sorted.set(index);
            }
        }
        return sorted;
    }
    
    /**
     * @return index of first set bit, or -1 if none set.
     */
    public static int first(BitVector v) {
        return first(v, 0);
    }
    
    /**
     * @return index of first set bit from starting index, or -1 if none set.
     */
    public static int first(BitVector v, int from) {
        Preconditions.checkArgument(from < v.size(), "From index must be less than size of vector.");
        BitVector subVector = BitUtils.newFromTo(v, from, v.size() - 1);
        if (subVector.cardinality() != 0) {
            for (int i = 0; i < subVector.size(); i++) {
                if (subVector.getQuick(i)) {
                    return i;
                }
            }
        }
        return -1;
    }
    
    public static List<Integer> assertedIndices(BitVector v, int from) {
        Preconditions.checkArgument(from < v.size(), "From index must be less than size of vector.");
        List<Integer> indices = Lists.newArrayList();
        for (int i = from; i < v.size(); i ++) {
            if (v.get(i)) {
                indices.add(i);
            }
        }
        return indices;
    }
    
    /**
     * @return v, BitVector of length with all bits true.
     */
    public static BitVector constant(int length) {
        BitVector v = new BitVector(length);
        v.not();
        return v;
    }
}

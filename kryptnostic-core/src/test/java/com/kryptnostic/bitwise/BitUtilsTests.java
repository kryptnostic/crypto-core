package com.kryptnostic.bitwise;

import java.util.List;
import java.util.Random;

import org.junit.Assert;
import org.junit.Test;

import cern.colt.bitvector.BitVector;

import com.kryptnostic.linear.BitUtils;

public class BitUtilsTests {
    private static final Random r = new Random(0);

    @Test
    public void testRandomVectorWithSpecificWeight() {
        int hammingWeight = r.nextInt(128);
        BitVector v = BitUtils.randomVector(128, hammingWeight);
        Assert.assertEquals(hammingWeight, v.cardinality());
    }

    @Test
    public void testRandomVectorReturnsCorrectLength() {
        int length = r.nextInt(128);
        BitVector v = BitUtils.randomVector(length);
        Assert.assertEquals(length, v.size());
    }

    @Test
    public void testExtend() {
        long[] values = { r.nextLong() };
        BitVector expected = new BitVector(values, 64);
        BitVector extended = BitUtils.extend(expected, 128);
        for (int i = 0; i < expected.size(); i++) {
            Assert.assertEquals(expected.get(i), extended.get(i));
        }
    }

    @Test
    public void testSortByMapping() {
        long[] values = { r.nextLong() };
        BitVector original = new BitVector(values, 64);
        int[] mapping = new int[original.size()];
        for (int i = 0; i < mapping.length; i++) {
            int index = ( i + 10 ) % mapping.length;
            mapping[i] = index;
        }

        BitVector sorted = BitUtils.extendAndOrder(original, mapping, original.size());

        for (int i = 0; i < original.size(); i++) {
            Assert.assertEquals(original.get(i), sorted.get(mapping[i]));
        }
    }
    
    @Test
    public void testFirst() {
        BitVector empty = new BitVector(64);
        Integer first = BitUtils.first(empty, 0);
        Assert.assertTrue(first == -1);
        
        for (int i = 0; i < 2000; i++) {
            BitVector random = BitUtils.randomVector(128);
            first = BitUtils.first(random, 0);
            Assert.assertTrue(random.get(first));
            if (first > 1) {
                Assert.assertTrue(BitUtils.newFromTo(random, 0, first-1).cardinality() == 0);
            } else if (first == 1) {
                Assert.assertFalse(random.get(0));
            }
        }
    }
    
    @Test
    public void testNewFromTo() {
        BitVector original = new BitVector(128);
        original.set(4);
        original.set(27);
        original.set(126);
        BitVector newVector = BitUtils.newFromTo(original, 0, 63);
        Assert.assertEquals(64, newVector.size());
        Assert.assertTrue(newVector.get(4));
        Assert.assertTrue(newVector.get(27));
        newVector.clear(4);
        newVector.clear(27);
        Assert.assertTrue(newVector.cardinality() == 0);
        
        newVector = BitUtils.newFromTo(original, 25, 127);
        Assert.assertEquals(103, newVector.size());
        Assert.assertTrue(newVector.get(2));
        Assert.assertTrue(newVector.get(101));
        newVector.clear(2);
        newVector.clear(101);
        Assert.assertTrue(newVector.cardinality() == 0);
    }
    
    @Test
    public void testAssertedIndices() {
        BitVector random = BitUtils.randomVector(128);
        List<Integer> asserted = BitUtils.assertedIndices(random, 0);
        for (Integer index : asserted) {
           Assert.assertTrue(random.get(index));
           random.clear(index);
        }
        Assert.assertTrue(random.cardinality() == 0);
    }
}

package com.kryptnostic.linear;

import java.util.Random;

import cern.colt.bitvector.BitVector;

public class BitUtils {
    private static final Random r = new Random( System.currentTimeMillis() );
    
    public static long parity( long l ) {
        l ^= l>>32;
        l ^= l>>16;
        l ^= l>>8;
        l ^= l>>4;
        l ^= l>>2;
        l ^= l>>1;
        return l & 1L;
    }
    
    public static BitVector randomBitVector( int length ) {
        //TODO: Optimize this
        BitVector v = new BitVector( length );
        for( int i = 0 ; i < length ; ++i ) {
            if( r.nextBoolean() ) {
                v.set( i );
            }
        }
        return v;
    }
 }

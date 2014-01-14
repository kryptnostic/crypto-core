package com.kryptnostic.multivariate;

import java.security.InvalidParameterException;
import java.util.Arrays;
import java.util.Random;
import java.util.Set;

import com.google.common.base.Preconditions;
import com.google.common.collect.Sets;

import cern.colt.bitvector.BitVector;

public class Monomial extends BitVector {
    private static final long serialVersionUID = -8751413919025034976L;

    public Monomial( int size ) {
        super( size );
    }
    
    public Monomial( long[] bits , int size ) {
        super( bits , size );
    }
    
    public boolean eval( BitVector input ) {
        if( size() == input.size() ) {
            BitVector check = copy();
            check.and( input );
            return check.equals( this ) ;
        } else {
            throw new InvalidParameterException("Number of terms in input doesn't not much number of terms in Monomial.");
        }
    }
    
    public Monomial product( Monomial monomial ) {
        Preconditions.checkArgument( this.size() == monomial.size() , "Cannot compute product due to polynomial ring mismatch.");
        Monomial result = new Monomial( monomial.size() );
        result.elements( Arrays.copyOf( this.elements() , this.elements().length ) , this.elements().length );
        result.and( monomial );
        return result;
    }
    
    public static Monomial randomMonomial( int size , int maxOrder ) {
        Random r = new Random( System.currentTimeMillis() );
        int order = r.nextInt( maxOrder - 1 ) + 1;
        Monomial monomial = new Monomial( size );
        
        Set<Integer> terms = Sets.newHashSet();
        while( terms.size() < order ) {
            terms.add( r.nextInt( size ) );
        }
        
        for( int term : terms ) {
            monomial.set( term );
        }
        
        return monomial;
    }
    
    public static Monomial constantMonomial( int size ) {
        return new Monomial( size );
    }
    
    public static Monomial linearMonomial( int size , int term ) {
        Monomial m = new Monomial( size );
        m.set( term );
        return m;
    }
}

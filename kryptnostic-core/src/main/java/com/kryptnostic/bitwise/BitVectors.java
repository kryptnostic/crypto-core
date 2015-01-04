package com.kryptnostic.bitwise;

import java.nio.ByteBuffer;
import java.nio.LongBuffer;
import java.security.SecureRandom;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.apache.commons.codec.binary.Base64;

import cern.colt.bitvector.BitVector;

import com.google.common.base.Function;
import com.google.common.base.Preconditions;
import com.google.common.collect.Iterables;
import com.google.common.collect.Sets;
import com.kryptnostic.linear.EnhancedBitMatrix;
import com.kryptnostic.multivariate.MultivariateUtils;

/**
 * Routines for manipulating bitvectors.
 * 
 * @author Matthew Tamayo-Rios &lt;matthew@kryptnostic.com&gt;
 * @author Sina Iman &lt;sina@kryptnostic.com&gt;
 * @author Nick Hewitt &lt;nick@kryptnostic.com&gt;
 */
public final class BitVectors {
    private static final Random r             = new SecureRandom();
    private static final int    INTEGER_BYTES = Integer.SIZE / Byte.SIZE;
    private static final Base64 codec         = new Base64();

    private BitVectors() {}

    private final static Function<BitVector, BitVector> cloner = new Function<BitVector, BitVector>() {
                                                                   @Override
                                                                   public BitVector apply( BitVector input ) {
                                                                       return input.copy();
                                                                   }
                                                               };

    /**
     * Returns an iterable that allow lazy evaluation of the cloning functions for efficient use of Guava collection
     * factory methods.
     * 
     * @param vectors
     * @return
     */
    public static Iterable<BitVector> cloneToIterable( BitVector... vectors ) {
        return Iterables.transform( Arrays.asList( vectors ), cloner );
    }

    /**
     * Performs a deep clone of Monomial class.
     * 
     * @param vectors to clone
     * @return An array of monomials generated by calling {@link BitVector#clone()} on each element in {@code monomials}
     */
    public static BitVector[] deepCloneBitvectorArray( BitVector... vectors ) {
        BitVector[] copies = new BitVector[ vectors.length ];
        for ( int i = 0; i < vectors.length; ++i ) {
            copies[ i ] = vectors[ i ].copy();
        }
        return copies;
    }

    /**
     * Turns a BitVector into a Base64 encoded string
     * 
     * @param input BitVector to be marshaled
     * @return Base64 encoded byte array in the following format: { bit_vector_size:int, bit_vector_bits[0]:long,
     *         bit_vector_bits[1]:long, ..., bit_vector_bits[n]:long }
     */
    public static String marshalBitvector( BitVector input ) {
        if ( input == null ) {
            return null;
        }
        long[] data = input.elements();
        byte[] target = new byte[ ( data.length << 3 ) + INTEGER_BYTES ];
        ByteBuffer buf = ByteBuffer.wrap( target );
        buf.putInt( input.size() );
        buf.asLongBuffer().put( data );
        return new String( codec.encode( target ) );
    }

    /**
     * Creates a BitVector from a Base64 encoded string
     * 
     * @param input Base64 encoded string of a byte array in the following format: { bit_vector_size:int,
     *            bit_vector_bits[0]:long, bit_vector_bits[1]:long, ..., bit_vector_bits[n]:long }
     * @return The unmarshaled BitVector
     */
    public static BitVector unmarshalBitvector( String input ) {
        if ( input == null ) {
            return null;
        }
        byte[] decoded = Base64.decodeBase64( input.getBytes() );
        ByteBuffer buf = ByteBuffer.wrap( decoded );
        int size = buf.getInt();
        LongBuffer longBuffer = buf.asLongBuffer();
        long[] longs = new long[ ( decoded.length - INTEGER_BYTES ) >>> 3 ];
        longBuffer.get( longs );
        return new BitVector( longs, size );
    }

    /**
     * Creates a bitvector of from an array of bytes. Bytes in excess of {@code bits} are truncated.
     * 
     * @param bits the total number of bits in the resultant bit vector. There must be enough byte supplied to create
     *            the bitvector, i.e {@code (bytes.length * 8) > bits}.
     * @param bytes an array of bytes used to create the bitvector
     * @return A bitvector of length {@code bits} where the i-th bit is {@code bytes[i/8] & (1<<(i %8))}.
     */
    public static BitVector fromBytes( int bits, byte[] bytes ) {
        ByteBuffer buf = ByteBuffer.wrap( bytes );
        BitVector result = new BitVector( bits );
        buf.asLongBuffer().get( result.elements() );
        return result;
    }

    /**
     * Concatenates two bitvectors.
     * 
     * @param lhs Lower order bitvector
     * @param rhs Higher order bitvector
     * @return Returns the bits from lhs as the lower order bits, followed by the bits from rhs as the higher order
     *         bits.
     */
    public static BitVector concatenate( BitVector lhs, BitVector rhs ) {
        BitVector result = new BitVector( lhs.size() + rhs.size() );
        result.replaceFromToWith( 0, lhs.size() - 1, lhs, 0 );
        result.replaceFromToWith( lhs.size(), lhs.size() + rhs.size() - 1, rhs, 0 );

        return result;
    }

    /**
     * Concatenates, in order, an array of BitVectors.
     * 
     * @param vectors The vectors to concatenate.
     * @return The vectors concatenated, in order, i.e {@code vectors[0],vectors[1],...,vectors[vectors.length -1]}}
     */
    public static BitVector concatenate( BitVector... vectors ) {
        return concatenate( Arrays.asList( vectors == null ? new BitVector[ 0 ] : vectors ) );
    }

    /**
     * Concatenates, in order, an array of BitVectors.
     * 
     * @param vectors The vectors to concatenate.
     * @return The vectors concatenated, in order, i.e {@code
     *         vectors.get(0),vectors.get(1),...,vectors.get(vectors.size()-1)}}
     */
    public static BitVector concatenate( List<BitVector> vectors ) {
        int baseIndex = 0;
        BitVector result = new BitVector( sumLengths( vectors ) );
        for ( BitVector v : vectors ) {
            result.replaceFromToWith( baseIndex, ( baseIndex += v.size() ) - 1, v, 0 );
        }
        return result;
    }

    /**
     * Extracts a bitvector from a matrix by concatenating its rows.
     * 
     * @param m The matrix m, whose rows will be concatenated. It must be a square n x n matrix.
     * @return The concatnated rows of {@code m}
     */
    public static BitVector fromSquareMatrix( EnhancedBitMatrix m ) {
        Preconditions.checkArgument( m.rows() == m.cols(), "Matrix must be square to derive bit vector from matrix." );
        return fromMatrix( m );
    }

    /**
     * Extracts a bitvector from a matrix by concatenating its rows.
     * 
     * @param m The matrix m, whose rows will be concatenated.
     * @return The concatnated rows of {@code m}
     */
    public static BitVector fromMatrix( EnhancedBitMatrix m ) {
        return concatenate( m.getRows() );
    }

    /**
     * Computes the sum of the lengths of a list of bitvectors.
     * 
     * @param vectors the bitvectors whose length will be totaled.
     * @return The sum of the length of each {@code BitVector} in {@code vectors}
     */
    public static int sumLengths( List<BitVector> vectors ) {
        int total = 0;
        for ( BitVector v : vectors ) {
            total += v.size();
        }
        return total;
    }

    /**
     * @param from , the index in the backing long array to start from
     * @param to , the last index in the backing long array to copy
     * @return BitVector
     */
    // TODO consider refactoring this so as not to reach into implementation.
    public static BitVector subVector( BitVector v, int from, int to ) {
        return new BitVector( Arrays.copyOfRange( v.elements(), from, to ), ( to - from ) << 6 );
    }

    public static BitVector randomVector( int length, int desiredHammingWeight ) {
        BitVector v = new BitVector( length );
        /*
         * In theory optimized popcnt instruction is going to be faster than bit twiddling to check individual bits.
         */
        while ( v.cardinality() < desiredHammingWeight ) {
            v.set( r.nextInt( length ) );
        }
        return v;
    }

    public static BitVector extend( BitVector v, int newSize ) {
        Preconditions.checkArgument( v.size() <= newSize, "New size must be greater than input vector size." );
        BitVector result = v.copy();
        result.setSize( newSize );
        return result;
    }

    public static BitVector randomVector( int length ) {
        return MultivariateUtils.randomVector( length );
    }

    /**
     * Given a mapping from old indices to new indices, creates a new BitVector ordered by this mapping.
     * 
     * @return BitVector
     */
    public static BitVector extendAndOrder( BitVector v, int[] mapping, int newSize ) {
        Preconditions.checkArgument(
                mapping.length == v.size(),
                "Must map exactly every bit in the vector to a new index." );
        Set<Integer> elements = Sets.newHashSet();
        BitVector sorted = new BitVector( newSize );
        for ( int i = 0; i < v.size(); i++ ) {
            int index = mapping[ i ];
            Preconditions.checkArgument( !elements.contains( index ), "Cannot map two variables to the same index." );
            elements.add( index );
            Preconditions.checkArgument( index < newSize, "Cannot map to index greater than new size." );
            if ( v.get( i ) ) {
                sorted.set( index );
            }
        }
        return sorted;
    }

    public static int getFirstSetBit( BitVector v ) {
        // TODO: Optimize
        for ( int i = 0; i < v.size(); ++i ) {
            if ( v.get( i ) ) {
                return i;
            }
        }
        return -1;
    }

    public static String asBitString( BitVector input ) {
        return asBitString( input, " " );
    }

    public static String asBitString( BitVector input, String separator ) {
        StringBuilder bitstring = new StringBuilder();
        int skipIndex = input.size() - 1;
        for ( int i = 0; i < input.size(); ++i ) {
            
            bitstring.append( input.get( i ) ? 1 : 0 );
            
            if( i != skipIndex ) {
                bitstring.append( separator );
            }
        }

        return bitstring.toString().trim();
    }
}

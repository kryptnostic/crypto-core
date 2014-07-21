package com.kryptnostic.crypto.fhe;

import java.util.Arrays;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.kryptnostic.crypto.PrivateKey;
import com.kryptnostic.linear.EnhancedBitMatrix;
import com.kryptnostic.linear.EnhancedBitMatrix.SingularMatrixException;
import com.kryptnostic.multivariate.PolynomialFunctionGF2;
import com.kryptnostic.multivariate.PolynomialFunctions;
import com.kryptnostic.multivariate.gf2.Monomial;
import com.kryptnostic.multivariate.gf2.PolynomialFunction;
import com.kryptnostic.multivariate.gf2.SimplePolynomialFunction;

public class HomomorphicFunctions {
    private static final Logger logger = LoggerFactory.getLogger( HomomorphicFunctions.class );
    public static SimplePolynomialFunction HomomorphicXor( int length , PrivateKey privateKey ) {
        return privateKey.computeHomomorphicFunction( PolynomialFunctions.XOR( length ) );
    }
    
    public static SimplePolynomialFunction HomomorphicAnd( int length , PrivateKey privateKey ) {
        return privateKey.computeHomomorphicFunction( PolynomialFunctions.AND( length ) );
    }
    
    public static SimplePolynomialFunction HomomorphicLsh( int length, PrivateKey privatekey ) {
        return privatekey.computeHomomorphicFunction( PolynomialFunctions.LSH( length , 1 ) );
    }
    
    public static SimplePolynomialFunction BinaryHomomorphicXor( int length , PrivateKey privateKey ) {
        return privateKey.computeBinaryHomomorphicFunction( PolynomialFunctions.BINARY_XOR( length ) );
    }
    
    public static SimplePolynomialFunction BinaryHomomorphicAnd( int length , PrivateKey privateKey ) {
        return privateKey.computeBinaryHomomorphicFunction( PolynomialFunctions.BINARY_AND( length ) );
    }

    public static SimplePolynomialFunction BinaryHomomorphicCarry( int length, PrivateKey privatekey ) {
        return privatekey.computeBinaryHomomorphicFunction( 
                PolynomialFunctions
                .LSH( length , 1 )
                .compose( PolynomialFunctions.BINARY_AND( length ) ) );
    }
    
    public static SimplePolynomialFunction HomomorphicHalfAdder( int length , PrivateKey privateKey ) {
        SimplePolynomialFunction xor = privateKey.computeBinaryHomomorphicFunction( PolynomialFunctions.BINARY_XOR( length ) );
        logger.info("Generated XOR portion of half adder.");
        SimplePolynomialFunction and = privateKey.computeBinaryHomomorphicFunction( PolynomialFunctions.BINARY_AND( length ) );
        logger.info("Generated AND portion of half adder.");
        SimplePolynomialFunction carry = 
                privateKey.computeBinaryHomomorphicFunction( 
                        PolynomialFunctions
                            .LSH( length , 1 )
                            .compose( PolynomialFunctions.BINARY_AND( length ) ) );
//        return privateKey.computeBinaryHomomorphicFunction( PolynomialFunctions.HALF_ADDER( 64 ) ) ; 
        logger.debug( "Generated carry portion of half adder" );
        return PolynomialFunctions.concatenate( xor , carry ); 
    }
    
    public static SimplePolynomialFunction DirectHomomorphicAnd( PrivateKey privateKey ) {
        SimplePolynomialFunction decryptor = privateKey.getDecryptor();
        Monomial [] monomials = decryptor.getMonomials();
        Monomial [] lhsMonomials = new Monomial[ monomials.length ];
        Monomial [] rhsMonomials = new Monomial[ monomials.length ];
        int inputLength = monomials[ 0 ].size() << 1;
        int outputLength = decryptor.getContributions()[ 0 ].size();
        
        for( int i = 0 ; i < monomials.length ; ++i ) {
            Monomial m = monomials[ i ];
            Monomial mLHS = new Monomial( Arrays.copyOf( m.elements() , m.elements().length << 1 ) , inputLength );
            Monomial mRHS = new Monomial( inputLength );
            long[] srcArray = m.elements();
            long[] destArray = mRHS.elements();
            for( int j = 0 ; j < srcArray.length ; ++j ) {
                destArray[ j + srcArray.length ] = srcArray[ j ];
            }
            lhsMonomials[ i ] = mLHS;
            rhsMonomials[ i ] = mRHS;
        }
        
        
        SimplePolynomialFunction X = new PolynomialFunctionGF2( inputLength , outputLength , lhsMonomials , decryptor.getContributions() );
        SimplePolynomialFunction Y = new PolynomialFunctionGF2( inputLength , outputLength , rhsMonomials , decryptor.getContributions() );
        logger.info("Generated functions for producting.");
        SimplePolynomialFunction XY = X.and( Y );
        logger.info("Computed product of decryption functons");
        
        return privateKey.encryptBinary( XY );
    }
    
    public static PolynomialFunction EfficientAnd( PrivateKey privateKey ) throws SingularMatrixException {
        EnhancedBitMatrix L  = privateKey.randomizedL(),
                          E1 = privateKey.getE1(),
                          E2 = privateKey.getE2(),
                          D  = privateKey.getD();
        int plaintextLength = E1.cols(),
            ciphertextLength = E1.rows();
        // TODO clean up declarations
        SimplePolynomialFunction X = PolynomialFunctions.lowerBinaryIdentity( ciphertextLength ); // matrix dimension mismatch
        SimplePolynomialFunction Y = PolynomialFunctions.upperBinaryIdentity( ciphertextLength );
        
        SimplePolynomialFunction DX = D.multiply( X );
        SimplePolynomialFunction DY = D.multiply( Y );
        SimplePolynomialFunction DXplusY = D.multiply( X.xor( Y ) );
        
        SimplePolynomialFunction F = privateKey.getF();
        SimplePolynomialFunction FofDX = F.compose( DX );
        SimplePolynomialFunction FofDXplusY = F.compose( DXplusY );
        
        EnhancedBitMatrix R = EnhancedBitMatrix.randomInvertibleMatrix( E1.rows() );
        SimplePolynomialFunction R1 = PolynomialFunctions.randomFunction( ciphertextLength << 1,  plaintextLength ) ,
                                 R2 = PolynomialFunctions.randomFunction( ciphertextLength << 1,  plaintextLength );
        SimplePolynomialFunction R1ofXY = R1.compose(X, Y);
        SimplePolynomialFunction R2ofXY = R2.compose(X, Y);
        
        SimplePolynomialFunction Lx = L.multiply( X );
        SimplePolynomialFunction Ly = L.multiply( Y );
        
        SimplePolynomialFunction V1 = 
                E1
                    .multiply( Lx.xor( R1ofXY ) )
                    .xor( E2.multiply( DXplusY.xor( R1ofXY ) ) );
        
        SimplePolynomialFunction V2 = 
                E1
                    .multiply( Ly.xor( R2ofXY ) )
                    .xor( E2.multiply( DXplusY.xor( R2ofXY ) ) );        
        
        SimplePolynomialFunction V3 = 
                E1
                    .multiply( R.multiply( FofDX.xor( R1ofXY ) ) )
                    .xor( E2.multiply( DXplusY.xor( R1ofXY ) ) );
        
        SimplePolynomialFunction V4 = 
                E1
                .multiply( R.multiply( FofDX.xor( R2ofXY ) ) )
                .xor( E2.multiply( DXplusY.xor( R2ofXY ) ) );
        
        
        SimplePolynomialFunction PLL =
                E1
                    .multiply( Lx.and( Ly ).xor( FofDXplusY ) ) 
                    .xor( E2.multiply( DXplusY ) );
        SimplePolynomialFunction PRL =
        		E1
        			.multiply( R.inverse().multiply( Lx ).and( Ly ).xor( FofDXplusY ) )
        			.xor( E2.multiply( DXplusY ) );
        SimplePolynomialFunction PLR = 
        		E1
        			.multiply( Lx.and( R.inverse().multiply( Ly ) ).xor( FofDXplusY ) )
        			.xor( E2.multiply( DXplusY ) );
        SimplePolynomialFunction PRR =
        		E1
	        		.multiply( R.inverse().multiply( Lx ).and( R.inverse().multiply( Ly ) ).xor( FofDXplusY ) ) 
	                .xor( E2.multiply( DXplusY ) );
        logger.info("Generated functions for producting.");
        
        SimplePolynomialFunction xor = PolynomialFunctions.XOR( ciphertextLength );
        SimplePolynomialFunction homomorphicXor = privateKey.computeHomomorphicFunction( xor );
        
        SimplePolynomialFunction PLLv1v2 = PLL.compose(V1, V2);
        SimplePolynomialFunction PRLv3v2 = PRL.compose(V3, V2);
        SimplePolynomialFunction PLRv1v4 = PLR.compose(V1, V4);
        SimplePolynomialFunction PRRv3v4 = PRR.compose(V3, V4);
        logger.info("Generated product parts.");
        
        return homomorphicXor.compose( 
        		homomorphicXor.compose( PLLv1v2, PRLv3v2 ), 
        		homomorphicXor.compose( PLRv1v4, PRRv3v4 ) );
    }
}

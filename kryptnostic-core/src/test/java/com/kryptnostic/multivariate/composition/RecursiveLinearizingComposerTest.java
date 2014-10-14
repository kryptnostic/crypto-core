package com.kryptnostic.multivariate.composition;

import java.util.concurrent.TimeUnit;

import org.junit.Assert;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import cern.colt.bitvector.BitVector;

import com.google.common.base.Stopwatch;
import com.kryptnostic.bitwise.BitVectors;
import com.kryptnostic.linear.EnhancedBitMatrix;
import com.kryptnostic.multivariate.BasePolynomialFunction;
import com.kryptnostic.multivariate.PolynomialFunctions;
import com.kryptnostic.multivariate.PolynomialFunctionsTests;
import com.kryptnostic.multivariate.gf2.SimplePolynomialFunction;

public class RecursiveLinearizingComposerTest {
    private static final Logger logger = LoggerFactory.getLogger( PolynomialFunctionsTests.class );
    
    @Test 
    public void testRecurisveLinearizingComposer() {
        BasePolynomialFunction f = (BasePolynomialFunction)PolynomialFunctions.denseRandomMultivariateQuadratic( 128 , 128);
        BasePolynomialFunction inner = (BasePolynomialFunction) EnhancedBitMatrix.randomMatrix( 128 , 256 ) .multiply(  PolynomialFunctions.identity( 256 ) );
        
        
        RecursiveLinearizingComposer composer = new RecursiveLinearizingComposer( f );
        Stopwatch watch = Stopwatch.createStarted();
        SimplePolynomialFunction composed = composer.compose( inner );
        logger.info( "Compose time: {} ms" , watch.elapsed( TimeUnit.MILLISECONDS ) );
        
        BitVector input =  BitVectors.randomVector( inner.getInputLength() );
        BitVector expected = f.apply( inner.apply( input ) );
        BitVector actual = composed.apply( input );
        
        Assert.assertEquals( expected , actual );
    }
}
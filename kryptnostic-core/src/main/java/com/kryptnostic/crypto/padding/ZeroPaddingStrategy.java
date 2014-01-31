package com.kryptnostic.crypto.padding;

import java.util.Arrays;

public class ZeroPaddingStrategy implements PaddingStrategy {
    public final int inputByteLength;
    public ZeroPaddingStrategy( int inputByteLength ) {
        this.inputByteLength = inputByteLength;
    }
    @Override
    public byte[] pad( byte[] unpadded ) {
//        int mask = unpadded.length & 7;
//        if( mask == 0 && (~mask)!=0 ) {
        int leftover = unpadded.length % inputByteLength;
        if( leftover != 0 ) {
            return Arrays.copyOf( unpadded , unpadded.length + inputByteLength - leftover );
        } else {
            return unpadded;
        }
    }
}

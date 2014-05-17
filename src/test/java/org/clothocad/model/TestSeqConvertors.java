/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.clothocad.model;

import org.clothocad.core.schema.SequenceConvertersTest;

/**
 *
 * @author prashantvaidyanathan
 */
public class TestSeqConvertors {
    public static void main(String[] args) {
        //SimpleSequenceTest.testFromNucSeq();
        SequenceConvertersTest x = new SequenceConvertersTest();
        x.testConvertNucSeqToSimpleSeq();
    }
    
}

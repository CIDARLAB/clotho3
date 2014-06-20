/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.clothocad.model;

import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.clothocad.core.schema.SequenceConvertersTest;

/**
 *
 * @author prashantvaidyanathan
 */
public class TestSeqConvertors {
    public static void main(String[] args) {
       
        try {
            SimpleSequenceTest x = new SimpleSequenceTest();
            x.testCanConvert();
            x.testConvertsTo();
            x.testConvert();
        } catch (IOException ex) {
            Logger.getLogger(TestSeqConvertors.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    
}

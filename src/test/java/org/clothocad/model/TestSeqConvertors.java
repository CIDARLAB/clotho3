/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.clothocad.model;

import java.io.IOException;
import java.lang.reflect.Method;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.clothocad.core.schema.ConverterMapper;
import org.clothocad.core.schema.RunScripts;

/**
 *
 * @author prashantvaidyanathan
 */
public class TestSeqConvertors {
    public static void main(String[] args) throws Exception {
       
        
        /*Class[] parameterTypes = new Class[1];
        parameterTypes[0] = String.class;
        Method method1 = TestSeqConvertors.class.getMethod("method1", parameterTypes);

        TestSeqConvertors demo = new TestSeqConvertors();
        demo.method2(demo, method1, "Hello World");*/
        testclasstester();
        
        try {
            SimpleSequenceTest x = new SimpleSequenceTest();
            x.testCanConvert();
            x.testConvertsTo();
            x.testConvert();
        } catch (IOException ex) {
            Logger.getLogger(TestSeqConvertors.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    
    public static void testclasstester() throws Exception
    {
        TestSeq tseq = new TestSeq("ATTAATTCT");
        ConverterMapper convmap = new ConverterMapper();
        convmap.testclass(tseq,"sequence");
    }
    
}

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.clothocad.core.schema;

import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.model.NucSeq;
import org.clothocad.model.SimpleSequence;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author prashantvaidyanathan
 */
public class SequenceConvertersTest {
    
    public SequenceConvertersTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    /**
     * Test of naiveConvertNucSeqToSimpleSeq method, of class SequenceConverters.
     */
    @Test
    public void testNaiveConvertNucSeqToSimpleSeq_NucSeq() {
        System.out.println("naiveConvertNucSeqToSimpleSeq");
        NucSeq nucObj = null;
        SequenceConverters instance = null;
        SimpleSequence expResult = null;
        SimpleSequence result = instance.naiveConvertNucSeqToSimpleSeq(nucObj);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of naiveConvertNucSeqToSimpleSeq method, of class SequenceConverters.
     */
    @Test
    public void testNaiveConvertNucSeqToSimpleSeq_SimpleSequence() {
        System.out.println("naiveConvertNucSeqToSimpleSeq");
        SimpleSequence simpleObj = null;
        SequenceConverters instance = null;
        NucSeq expResult = null;
        NucSeq result = instance.naiveConvertNucSeqToSimpleSeq(simpleObj);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of convertNucSeqToSimpleSeq method, of class SequenceConverters.
     */
    @Test
    public void testConvertNucSeqToSimpleSeq() {
        try {
            System.out.println("convertNucSeqToSimpleSeq");
            NucSeq nucObj = new NucSeq("atttccccggggaccccccctttt");
            SimpleSequence expResult = null;
            SimpleSequence result = SequenceConverters.convertNucSeqToSimpleSeq(nucObj);
            assertEquals(expResult, result);
            // TODO review the generated test code and remove the default call to fail.
            fail("The test case is a prototype.");
        } catch (IllegalArgumentException ex) {
            Logger.getLogger(SequenceConvertersTest.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            Logger.getLogger(SequenceConvertersTest.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    /**
     * Test of guardedConvert method, of class SequenceConverters.
     */
    @Test
    public void testGuardedConvert_Map_Schema() {
        System.out.println("guardedConvert");
        Map data = null;
        Schema type = null;
        SequenceConverters instance = null;
        ObjBase expResult = null;
        ObjBase result = instance.guardedConvert(data, type);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of guardedConvert method, of class SequenceConverters.
     */
    @Test
    public void testGuardedConvert_Map_String() {
        System.out.println("guardedConvert");
        Map data = null;
        String schemaName = "";
        SequenceConverters instance = null;
        ObjBase expResult = null;
        ObjBase result = instance.guardedConvert(data, schemaName);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }
    
}

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.clothocad.model;

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
public class SimpleSequenceTest {
    
    public SimpleSequenceTest() {
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
     * Test of fromNucSeq method, of class SimpleSequence.
     */
  
    /**
     * Test of getSequence method, of class SimpleSequence.
     */
    @Test
    public void testGetSequence() {
        System.out.println("getSequence");
        SimpleSequence instance = null;
        String expResult = "";
        String result = instance.getSequence();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of setSequence method, of class SimpleSequence.
     */
    @Test
    public void testSetSequence() {
        System.out.println("setSequence");
        String sequence = "";
        SimpleSequence instance = null;
        instance.setSequence(sequence);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }
    
}

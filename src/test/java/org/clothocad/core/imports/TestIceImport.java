/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.imports;
import java.util.ArrayList;
import java.util.List;
import org.clothocad.core.layers.communication.imports.IceImporter;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.utils.TestUtils;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author spaige
 */
public class TestIceImport {
    
    private static Persistor p = new TestUtils().getA(Persistor.class);
    public TestIceImport() {
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

    @Test
    public void testIceImport(){
        List<Integer> ints = new ArrayList<>();
        ints.add(255);
        ints.add(415);
        
        IceImporter.directImport(p, ints, 0);
        
    }
    
    @Test 
    public void testIceFaster(){
        List<Integer> ints = new ArrayList<>();
        ints.add(255);
        ints.add(415);
        
        IceImporter.importData(ints);
    }    
}
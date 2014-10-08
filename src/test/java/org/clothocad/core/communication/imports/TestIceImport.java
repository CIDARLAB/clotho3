/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.communication.imports;
import java.util.ArrayList;
import java.util.List;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.util.TestUtils;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

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

    
    //Putting these on ice (ha ha), the ICE api has changed.
    
    @Ignore @Test
    public void testIceImport(){
        List<Integer> ints = new ArrayList<>();
        ints.add(255);
        ints.add(415);
        
        IceImporter.directImport(p, ints, 0);
        
    }
    
    @Ignore @Test 
    public void testIceFaster(){
        List<Integer> ints = new ArrayList<>();
        ints.add(255);
        ints.add(415);
        
        IceImporter.importData(ints);
    }    
}
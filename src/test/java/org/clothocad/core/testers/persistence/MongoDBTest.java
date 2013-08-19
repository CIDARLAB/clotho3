
package org.clothocad.core.testers.persistence;

import org.clothocad.model.FreeForm;
import com.google.inject.Guice;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertSame;

import java.io.IOException;
import java.util.Random;

import org.bson.types.ObjectId;
import org.clothocad.core.datums.ObjBase;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import com.mongodb.BasicDBObject;
import java.net.UnknownHostException;
import java.util.Collection;
import org.clothocad.core.persistence.ClothoConnection;
import org.clothocad.core.persistence.mongodb.MongoDBConnection;
import org.clothocad.core.utils.TestUtils;
import org.clothocad.model.Feature;
import org.clothocad.model.Institution;
import org.clothocad.model.Lab;
import org.clothocad.model.Part;
import org.clothocad.model.Person;
import org.junit.After;

public class MongoDBTest {

    static private ClothoConnection conn;
    
    @BeforeClass
    public static void setUpClass() throws UnknownHostException {
        conn = TestUtils.getDefaultTestInjector().getInstance(MongoDBConnection.class);
    }
    
    @Before
    public void setUp() {
        conn.deleteAll();
    }
    
    @After
    public void tearDown() {}
    
    public static void saveAndGet(ObjBase o){
        conn.save(o);
        ObjectId id = o.getUUID();
        Class c = o.getClass();
        o = null;
        conn.get(c, id);
    }
    
    @Test
    public void testSimpleObject() throws IOException{
        
        //simple data class
        Institution i = new Institution("Test institution", "Townsville", "Massachusetts", "United States of America");
        conn.save(i);
        ObjectId id = i.getUUID();
        i = null;
        BasicDBObject query = new BasicDBObject("name", "Test institution");
        i =  conn.get(Institution.class, id);
        assertEquals("Townsville",i.getCity());
        assertEquals("Massachusetts",i.getState());
        assertEquals("United States of America",i.getCountry());
    }
    
    @Test 
    public void testCompositeObjectThroughDBConnection() throws IOException{
        Institution i = new Institution("Test institution", "Townsville", "Massachusetts", "United States of America");
        Lab lab = new Lab(i, null, "Test Lab", "College of Testing", "8 West Testerfield");
        
        conn.save(i);
        conn.save(lab);
        ObjectId id = lab.getUUID();
        
        lab = null;
        
        lab = conn.get(Lab.class, id);
        
        assertEquals(i.getUUID(),lab.getInstitution().getUUID());
    }
    @Test
    public void testIdAssociation(){
        Institution i = new Institution("Test institution", "Townsville", "Massachusetts", "United States of America");
        assertNull(i.getUUID());
        conn.save(i);
        assertNotNull(i.getUUID());
    }
    
    @Test 
    public void testGetById(){
        Institution i = new Institution("Test institution", "Townsville", "Massachusetts", "United States of America");
        ObjectId id = ObjectId.get();
        i.setUUID(id);
        conn.save(i);
        Institution j = conn.get(Institution.class, id);
        assertEquals(i.getUUID(), j.getUUID());
        assertEquals(i.getCity(), j.getCity());
        assertEquals(i.getCountry(), j.getCountry());
        assertEquals(i.getState(), j.getState());
        
    }
    
    
    /*@Test
    public void testConflict(){
        //figure out what to do when overwriting more recent changes
    }*/
    
    
    //TODO
    public void testCreateFromJSON() {
        
    }
     
    //TODO
    public void testSimpleSuperAndSubClass(){
 
    }
}

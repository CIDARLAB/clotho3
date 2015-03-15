
package org.clothocad.core.persistence;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertFalse;

import java.io.IOException;

import org.clothocad.core.datums.ObjBase;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import com.mongodb.BasicDBObject;
import java.net.UnknownHostException;
import java.util.Map;
import static org.clothocad.core.ReservedFieldNames.*;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.jongo.JongoConnection;
import org.clothocad.core.util.TestUtils;
import org.clothocad.model.Institution;
import org.clothocad.model.Lab;
import org.junit.After;

public class MongoDBTest {

    static private ClothoConnection conn;
    
    @BeforeClass
    public static void setUpClass() throws UnknownHostException {
        conn = TestUtils.getDefaultTestInjector().getInstance(JongoConnection.class);
    }
    
    @Before
    public void setUp() {
        conn.deleteAll();
    }
    
    @After
    public void tearDown() {}
    
    public static void saveAndGet(ObjBase o){
        conn.save(o);
        ObjectId id = o.getId();
        Class c = o.getClass();
        o = null;
        conn.get(c, id);
    }
    
    @Test
    public void testSimpleObject() throws IOException{
        
        //simple data class
        Institution i = new Institution("Test institution", "Townsville", "Massachusetts", "United States of America");
        conn.save(i);
        ObjectId id = i.getId();
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
        ObjectId id = lab.getId();
        
        lab = null;
        
        lab = conn.get(Lab.class, id);
        
        assertEquals(i.getId(),lab.getInstitution().getId());
    }
    @Test
    //id is non-null after saving
    public void testIdAssociation(){
        Institution i = new Institution("Test institution", "Townsville", "Massachusetts", "United States of America");
        assertNull(i.getId());
        conn.save(i);
        assertNotNull(i.getId());
    }
    
    @Test 
    public void testGetById(){
        Institution i = new Institution("Test institution", "Townsville", "Massachusetts", "United States of America");
        ObjectId id = new ObjectId();
        i.setId(id);
        conn.save(i);
        Institution j = conn.get(Institution.class, id);
        assertEquals(i.getId(), j.getId());
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
    
    //Test that returned maps have  field instead of "_id"
    @Test
    public void testIdFieldName(){
        Institution i = new Institution("Test institution", "Townsville", "Massachusetts", "United States of America");
        ObjectId id = new ObjectId();
        i.setId(id);
        conn.save(i);
        Map<String,Object> map = conn.getAsBSON(id);
        assertTrue(map.containsKey(ID));
        assertFalse(map.containsKey("_id"));
    }
}

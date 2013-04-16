
package org.clothocad.core.testers.persistence;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.Random;

import org.bson.types.ObjectId;
import org.clothocad.core.aspects.Persistor;
import org.clothocad.core.datums.ObjBase;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import com.mongodb.BasicDBObject;
import com.mongodb.MongoClient;
import java.net.UnknownHostException;
import org.clothocad.core.layers.persistence.ClothoConnection;
import org.clothocad.core.layers.persistence.mongodb.MongoDBConnection;
import org.clothocad.model.Feature;
import org.clothocad.model.Institution;
import org.clothocad.model.Lab;
import org.clothocad.model.Part;
import org.clothocad.model.Person;
import org.junit.After;
import org.junit.AfterClass;

public class MongoDBTester {
    
    public MongoDBTester() {
    }
    
    static private MongoClient mongo;
    static private ClothoConnection conn;
    
    @BeforeClass
    public static void setUpClass() throws UnknownHostException {
        
        mongo = new MongoClient();
        conn = new MongoDBConnection();
    }
    
    @Before
    public void setUp() {
        mongo.dropDatabase("clotho");

    }
    
    @After
    public void tearDown() {

    }
    // TODO: tests that only test MongoDBConnection
    
    public static void saveAndGet(ObjBase o){
        o.save();
        ObjectId id = o.getUUID();
        Class c = o.getClass();
        o = null;
        conn.get(c, id);
    }
    
    @Test
    public void testSimpleObject() throws IOException{
        
        //simple data class
        Institution i = new Institution("Test institution", "Townsville", "Massachusetts", "United States of America");
        boolean status = i.save();
        assertTrue(status);
        ObjectId id = i.getUUID();
        i = null;
        BasicDBObject query = new BasicDBObject("name", "Test institution");
        i =  conn.get(Institution.class, id);
        assertEquals("Townsville",i.getCity());
        assertEquals("Massachusetts",i.getState());
        assertEquals("United States of America",i.getCountry());
    }
    
    @Test
    public void testSharedObject(){
        
        Person testPerson = new Person("Test Person", null, null);
        
        //class w/ composition
        Part part1 = Part.generateBasic("test part", "This part is a test", "ATCG", new FreeForm(), testPerson);
        Part part2 = Part.generateBasic("different test part", "This part is another test", "TCAG", new FreeForm(), testPerson);
        boolean status = part1.save();
        assertTrue(status);
        
        ObjectId id1 = part1.getUUID();
        
        //can now find testPerson in DB
        assertNotNull(testPerson.getUUID());
        

        testPerson.setDisplayName("Different Name");
        status = part2.save();
        assertTrue(status);
        
        
        part1 = conn.get(Part.class, id1);
        
        //changes in one composite cause changes in the other
        String name = part1.getAuthor().getDisplayName();
        assertEquals("Different Name",name);
    }
    
    @Test
    public void testSuperAndSubClass(){
        Part p = Part.generateBasic("test part", "This part is a test", "ATCG", new FreeForm(), null);
        p.save();
        
        ObjectId id = p.getUUID();
        ExtendedPart ep = conn.get(ExtendedPart.class, id);
        
        assertEquals(p.getUUID(), ep.getUUID());
        assertEquals("test part",ep.getName());
        assertEquals("This part is a test",ep.getShortDescription());
        
        ep.setAdditionalParameters("test params");
        
        ep.save();
        ep = null;
        p = conn.get(Part.class, id);
        p.setName("renamed part");
        p.save();
        
        ep = conn.get(ExtendedPart.class, id);
        assertEquals("test params",ep.getAdditionalParameters() );
        assertEquals("renamed part",ep.getName());
        
    }
    
    
    @Test 
    public void testCompositeObject() throws IOException{
        Institution i = new Institution("Test institution", "Townsville", "Massachusetts", "United States of America");
        Lab lab = new Lab(i, null, "Test Lab", "College of Testing", "8 West Testerfield");
        saveAndGet(lab);
        
    }
    @Test 
    public void testCompositeObjectThroughDBConnection() throws IOException{
        Institution i = new Institution("Test institution", "Townsville", "Massachusetts", "United States of America");
        Lab lab = new Lab(i, null, "Test Lab", "College of Testing", "8 West Testerfield");
        conn.save(i);
        conn.save(lab);
    }
    @Test
    public void testIdAssociation(){
        Institution i = new Institution("Test institution", "Townsville", "Massachusetts", "United States of America");
        assertNull(i.getUUID());
        i.save();
        assertNotNull(i.getUUID());
    }
    
    @Test 
    public void testGetById(){
        Institution i = new Institution("Test institution", "Townsville", "Massachusetts", "United States of America");
        ObjectId id = ObjectId.get();
        i.setUUID(id);
        i.save();
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
    
    @Test
    public void testCircular(){
        Institution i = new Institution("Test institution", "Townsville", "Massachusetts", "United States of America");
        Lab lab = new Lab(i, null, "Test Lab", "College of Testing", "8 West Testerfield");
        Person testPerson = new Person("Test Person", lab, null);
        lab.setPI(testPerson);
        
        saveAndGet(testPerson);
        
        assertEquals("Test Person", testPerson.getDisplayName());
        assertEquals("Test institution", testPerson.getLab().getInstitution().getName());
        assertSame(testPerson, testPerson.getLab().getPI());
    }
    
    @Test
    public void testCreateGFPInstance(){
        Feature gfp = Feature.generateFeature("GFPuv", "ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA", 
                    null, true);
        saveAndGet(gfp);
        
        
    }
    
    @Test
    public void testGetAllParts() {
    	int N = 100;

	   	// first, we create N parts
	   	for(int i=1; i<=N; i++) {
            (Part.generateBasic("part-"+i, "This is test part "+i, randomSequence(i), new FreeForm(), null))
            	.save();
    	}
    	
    	// then, retrieve all parts
    	Part[] arr = conn.get(Part.class);
    	assertEquals(arr.length, N); 
    }
    
    @Test
    public void testGetAllPartSequences() {    	
    	// how can I get all parts whose sequence is circular ??
    }
    
    private String randomSequence(int length) {
		Random randomGenerator = new Random();
		
		StringBuilder sb = new StringBuilder();
		for(int i=1; i<=10; i++) {
			int r = randomGenerator.nextInt(4);
			if(r == 0) {
				sb.append("A");
			} else if(r == 1) {
				sb.append("T");
			} else if(r == 2) {
				sb.append("C");
			} else if(r == 3) {
				sb.append("G");
			}
		}
		return sb.toString();
    } 
}

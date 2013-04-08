
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

import org.clothocad.model.Feature;
import org.clothocad.model.Institution;
import org.clothocad.model.Lab;
import org.clothocad.model.Part;
import org.clothocad.model.Person;

public class MongoDBTest {
    
    public static void saveAndGet(ObjBase o){
        o.save();
        ObjectId id = o.getUUID();
        Class c = o.getClass();
        o = null;
        Persistor.get().get(c, id);
    }
    
    @BeforeClass
    public static void clearDB() {    	
    	// clear the DB content before running the unit tests...
    	// so we can look into the DB content after the unit tests
    	Persistor.get().clearDB();
    }
    
    @Test
    public void testSimpleObject() throws IOException{
        
        //simple data class
        Institution i = new Institution("Test institution", "Townsville", "Massachusetts", "United States of America");
        boolean status = i.save();
        assertTrue(status);
        ObjectId id = i.getUUID();
        i = null;
        
        //BasicDBObject query = new BasicDBObject("name", "Test institution");  // ???
        
        i = (Institution)Persistor.get().get(Institution.class, id);
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
        
        part1 = (Part)Persistor.get().get(Part.class, id1);
        
        //changes in one composite cause changes in the other
        String name = part1.getAuthor().getDisplayName();
        assertEquals("Different Name",name);
    }
    
    @Test
    public void testSuperAndSubClass(){
        Part p = Part.generateBasic("test part", "This part is a test", "ATCG", new FreeForm(), null);
        p.save();
        
        ObjectId id = p.getUUID();
        ExtendedPart ep = (ExtendedPart)Persistor.get().get(ExtendedPart.class, id);
        
        assertEquals(p.getUUID(), ep.getUUID());
        assertEquals("test part",ep.getName());
        assertEquals("This part is a test",ep.getShortDescription());
        
        ep.setAdditionalParameters("test params");
        
        ep.save();
        ep = null;
        
        p = (Part)Persistor.get().get(Part.class, id);
        p.setName("renamed part");
        p.save();
        
        ep = (ExtendedPart)Persistor.get().get(ExtendedPart.class, id);
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
        i.save();

        Lab lab = new Lab(i, null, "Test Lab", "College of Testing", "8 West Testerfield");
        lab.save();
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
        Institution j = (Institution)Persistor.get().get(Institution.class, id);
        assertEquals(i.getUUID(), j.getUUID());
        assertEquals(i.getCity(), j.getCity());
        assertEquals(i.getCountry(), j.getCountry());
        assertEquals(i.getState(), j.getState());
        
    }
    
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
        Feature gfp = Feature.generateFeature("GFPuv", 
        		"ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA", 
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
    	Part[] arr = Persistor.get().get(Part.class);
    	assertEquals(arr.length, N+3);  // +3 since 3 parts have been inserted by the earlier tests
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

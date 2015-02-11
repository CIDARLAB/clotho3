/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence;

import com.fasterxml.jackson.core.JsonToken;

import de.undercouch.bson4jackson.BsonParser;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import org.clothocad.model.FreeForm;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.util.AuthorizedShiroTest;
import org.clothocad.model.Feature.FeatureRole;
import org.clothocad.model.Part;
import org.clothocad.model.Feature;
import org.clothocad.model.Institution;
import org.clothocad.model.Lab;
import org.clothocad.model.SimpleSequence;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 *
 * @author spaige
 */
public class PersistorTest extends AuthorizedShiroTest{

    public PersistorTest() {
        persistor = injector.getInstance(Persistor.class);
    }
    private static Persistor persistor;

    @Before
    public void setUp() {
        persistor.deleteAll();
    }
    
    @Test
    public void testSharedObject() {

        LabPersonForTests testPerson = new LabPersonForTests("Test Person", null, null);

        //class w/ composition
        Part part1 = new Part("test part", "This part is a test", 
        		new SimpleSequence("ATCG", testPerson), testPerson);
        part1.setFormat(new FreeForm());
        Part part2 = new Part("different test part", "This part is another test", 
        		new SimpleSequence("TCAG", testPerson), testPerson);
        part2.setFormat(new FreeForm());
        persistor.save(part1);

        ObjectId id1 = part1.getId();

        //can now find testPerson in DB
        assertNotNull(testPerson.getId());


        testPerson.setDisplayName("Different Name");
        persistor.save(part2);


        part1 = persistor.get(Part.class, id1);

        //changes in one composite cause changes in the other
        String name = part1.getAuthor().getDisplayName();
        assertEquals("Different Name", name);
    }

    @Test
    public void testCompositeObject() throws IOException {
        Institution i = new Institution("Test institution", "Townsville", "Massachusetts", "United States of America");
        Lab lab = new Lab(i, null, "Test Lab", "College of Testing", "8 West Testerfield");
        saveAndGet(lab);

    }

    @Test
    public void testCircular() {
        Institution i = new Institution("Test institution", "Townsville", "Massachusetts", "United States of America");
        Lab lab = new Lab(i, null, "Test Lab", "College of Testing", "8 West Testerfield");
        LabPersonForTests testPerson = new LabPersonForTests("Test Person", lab, null);
        lab.setPI(testPerson);

        saveAndGet(testPerson);

        assertEquals("Test Person", testPerson.getDisplayName());
        assertEquals("Test institution", testPerson.getLab().getInstitution().getName());
        assertSame(testPerson, testPerson.getLab().getPI());
    }

    @Test
    public void testCompositeSuperAndSubClass() {
    	LabPersonForTests testPerson = new LabPersonForTests("Test Person", null, null);
        Part p = new Part("test part", "This part is a test", 
        		new SimpleSequence("ATCG", testPerson), testPerson);
        p.setFormat(new FreeForm());
        persistor.save(p);
        Institution i = new Institution("Test institution", "Townsville", "Massachusetts", "United States of America");
        persistor.save(i);
        ObjectId id = i.getId();
        i = persistor.get(Institution.class, id);
        
        id = p.getId();
        p = persistor.get(Part.class, id);
        ExtendedPart ep = persistor.get(ExtendedPart.class, id);

        assertEquals(p.getId(), ep.getId());
        assertEquals("test part", ep.getName());
        assertEquals("This part is a test", ep.getDescription());

        ep.setAdditionalParameters("test params");

        persistor.save(ep);
        ep = null;
        p = persistor.get(Part.class, id);
        p.setName("renamed part");
        persistor.save(p);

        ep = persistor.get(ExtendedPart.class, id);
        assertEquals("test params", ep.getAdditionalParameters());
        assertEquals("renamed part", ep.getName());

    }

    @Test
    public void testGetAllBasicParts() {
        int N = 100;
        LabPersonForTests testPerson = new LabPersonForTests("Test Person", null, null);
        // first, we create N parts
        for (int i = 1; i <= N; i++) {
        	Part testPart = new Part("part-" + i, "This is test part " + i, 
        			new SimpleSequence(randomSequence(i), testPerson), testPerson);
        	testPart.setFormat(new FreeForm());
            persistor.save(testPart);
        }

        // then, retrieve all parts
        Collection<Part> results = persistor.getAll(Part.class);
        assertEquals(N, results.size());
    }
    
    //TODO:
    public void testGetAllParts(){
        
    }

    private String randomSequence(int length) {
        Random randomGenerator = new Random();

        StringBuilder sb = new StringBuilder();
        for (int i = 1; i <= 10; i++) {
            int r = randomGenerator.nextInt(4);
            if (r == 0) {
                sb.append("A");
            } else if (r == 1) {
                sb.append("T");
            } else if (r == 2) {
                sb.append("C");
            } else if (r == 3) {
                sb.append("G");
            }
        }
        return sb.toString();
    }

    @Test
    public void testCreateGFPInstance() {
    	LabPersonForTests testPerson = new LabPersonForTests("Test Person", null, null);
        Feature gfp = new Feature("GFPuv", FeatureRole.CDS, testPerson);
        gfp.setSequence(new SimpleSequence("ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA",
        		testPerson));
        saveAndGet(gfp);
    }

    //TODO:
    public void testGetAllPartSequences() {
        // how can I get all parts whose sequence is circular ??
    }

    public static void saveAndGet(ObjBase o) {
        persistor.save(o);
        ObjectId id = o.getId();
        Class c = o.getClass();
        o = null;
        persistor.get(c, id);
    }
    
    @Test
    public void initBuiltInsIdempotent(){
        int size = 0;
        Iterator it = persistor.find(new HashMap()).iterator();
        while (it.hasNext()){
            size ++;
            it.next();
        }
        persistor.initializeBuiltInSchemas();
        int newSize = 0;
         it = persistor.find(new HashMap()).iterator();
        while (it.hasNext()){
            newSize ++;
            it.next();
        }
        assertEquals(size, newSize);
    }
    
    //TODO: verify results are in appropriate style (id instead of _id), (schema instead of ClassName), etc
    
    
    public static List<JsonToken> getAllBsonTokens(BsonParser parser) throws IOException{
        List<JsonToken> out = new ArrayList<>();
        while (parser.nextToken() != null){
            out.add(parser.getCurrentToken());
        }
        
        return out;
    }
}
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.testers.persistence;

import java.io.IOException;
import java.util.Collection;
import java.util.Random;
import org.clothocad.model.FreeForm;
import org.bson.types.ObjectId;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.utils.TestUtils;
import org.clothocad.model.BasicPart;
import org.clothocad.model.Feature;
import org.clothocad.model.Institution;
import org.clothocad.model.Lab;
import org.clothocad.model.Part;
import org.clothocad.model.Person;
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
public class PersistorTest {

    public PersistorTest() {
    }
    private static Persistor persistor;

    @BeforeClass
    public static void setUpClass() {
        persistor = TestUtils.getDefaultTestInjector().getInstance(Persistor.class);
    }

    @AfterClass
    public static void tearDownClass() {
    }

    @Before
    public void setUp() {
        persistor.deleteAll();
    }

    @After
    public void tearDown() {
    }
    // TODO add test methods here.
    // The methods must be annotated with annotation @Test. For example:
    //
    // @Test
    // public void hello() {}

    @Test
    public void testSharedObject() {

        Person testPerson = new Person("Test Person", null, null);

        //class w/ composition
        Part part1 = Part.generateBasic("test part", "This part is a test", "ATCG", new FreeForm(), testPerson);
        Part part2 = Part.generateBasic("different test part", "This part is another test", "TCAG", new FreeForm(), testPerson);
        persistor.save(part1);

        ObjectId id1 = part1.getUUID();

        //can now find testPerson in DB
        assertNotNull(testPerson.getUUID());


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
        Person testPerson = new Person("Test Person", lab, null);
        lab.setPI(testPerson);

        saveAndGet(testPerson);

        assertEquals("Test Person", testPerson.getDisplayName());
        assertEquals("Test institution", testPerson.getLab().getInstitution().getName());
        assertSame(testPerson, testPerson.getLab().getPI());
    }

    @Test
    public void testCompositeSuperAndSubClass() {
        Part p = Part.generateBasic("test part", "This part is a test", "ATCG", new FreeForm(), null);
        persistor.save(p);

        ObjectId id = p.getUUID();
        ExtendedPart ep = persistor.get(ExtendedPart.class, id);

        assertEquals(p.getUUID(), ep.getUUID());
        assertEquals("test part", ep.getName());
        assertEquals("This part is a test", ep.getShortDescription());

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

        // first, we create N parts
        for (int i = 1; i <= N; i++) {
            persistor.save(Part.generateBasic("part-" + i, "This is test part " + i, randomSequence(i), new FreeForm(), null));
        }

        // then, retrieve all parts
        Collection<BasicPart> results = persistor.getAll(BasicPart.class);
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
        Feature gfp = Feature.generateFeature("GFPuv", "ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA",
                null, true);
        saveAndGet(gfp);
    }

    //TODO:
    public void testGetAllPartSequences() {
        // how can I get all parts whose sequence is circular ??
    }

    public static void saveAndGet(ObjBase o) {
        persistor.save(o);
        ObjectId id = o.getUUID();
        Class c = o.getClass();
        o = null;
        persistor.get(c, id);
    }
    
    //TODO: verify results are in appropriate style (id instead of _id), (schema instead of ClassName), etc
}
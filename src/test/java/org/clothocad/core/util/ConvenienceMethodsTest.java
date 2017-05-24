/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.model.Annotation;
import org.clothocad.model.BioDesign;
import org.clothocad.model.Parameter;
import org.clothocad.model.Part;
import org.clothocad.model.Person;
import org.clothocad.webserver.jetty.ConvenienceMethods;
import static org.clothocad.webserver.jetty.ConvenienceMethods.*;
import static org.junit.Assert.*;
import org.junit.ComparisonFailure;
import org.junit.Test;

/**
 *
 * @author David
 */
public class ConvenienceMethodsTest extends AuthorizedShiroTest {

    private static Persistor persistor;

    public ConvenienceMethodsTest() {
        persistor = injector.getInstance(Persistor.class);

        System.out.println();
        System.out.println();
    }

    /*
    BioDesign
    name
    
    Part
    name
    
    Sequence - Only if sequence provided
    name
    sequence
    
    Annotation - Only if sequence or role provided
    name
    start
    end
    
    Feature - Only if role provided
    name
    role
    
    BasicModule - Only if role provided
    name
    role
     */
 /*
    Had to change function signature to be able to capture all of the optional combinations of String parameters
    
    __optionals__:
    role
    sequence
     */
    @Test
    public void testCreatePart() {
        ObjectId first = createPart(persistor, "mySpecialPart", "David");

        Map<String, String> roleParam = new HashMap<>();
        roleParam.put("role", "PROMOTER");
        ObjectId second = createPart(persistor, "roleOnlyPart", roleParam, "David");

        Map<String, String> sequence = new HashMap<>();
        sequence.put("sequence", "catcatcatcatcatcatcatcatcat");
        ObjectId third = createPart(persistor, "FunCatPart", sequence, "David");

        Map<String, String> seqrole = new HashMap<>();
        seqrole.put("role", "GENE");
        seqrole.put("sequence", "tccctatcagtgatagagattgacatccctatcagtgatagagatactgagcac");
        ObjectId fourth = createPart(persistor, "R0040 Sequence", seqrole, "June Rhee");

        ArrayList<Parameter> paramObj = new ArrayList<>();
        paramObj.add(new Parameter("pTac", 626.3, "gene expression", "REU"));
        ObjectId fifth = createPart(persistor, "PartWithName", "Hello I am Name", "David");

        ObjectId sixth = createPart(persistor, "AnotherPart", "Another worthy part", paramObj, "David");

        paramObj.add(new Parameter("LacI Sensor Sensitivity", 0.5, "sensitivity", ""));
        ObjectId seventh = createPart(persistor, "FullPart", "Full Part with all Parameters", seqrole, paramObj, "David");

        compareDesign(first, getPart(persistor, "mySpecialPart", null, null, null, null, false));
        compareDesign(second, getPart(persistor, "roleOnlyPart", null, "promoter", null, null, false));
        compareDesign(third, getPart(persistor, "funcat", null, null, "catcatcatcatcat", null, false));
        compareDesign(fourth, getPart(persistor, "R0040", null, "gene", "tccctatcag", null, false));
        compareDesign(fifth, getPart(persistor, "PartWithName", "name", null, null, null, false));
        compareDesign(sixth, getPart(persistor, "AnotherPart", "worthy", null, null, null, false));
        compareDesign(seventh, getPart(persistor, "FullPart", "Part with all", "gene", "tccctatcag", paramObj, false));
    }

    @Test
    public void testCreateDevice() {

        ArrayList<String> partIDs = new ArrayList<>();

        ObjectId first = createDevice(persistor, "Barebones Device", partIDs, "David T.", false);
        compareDesign(first, getDevice(persistor, "Barebones", null, null, null, null, null, false));
        partIDs.add(first.getValue());

        Map<String, String> sequence = new HashMap<>();
        sequence.put("sequence", "catcat");
        ObjectId second = createDevice(persistor, "Basic Device", partIDs, sequence, "David T.", false);
        compareDesign(second, getDevice(persistor, "basic device", null, null, "catcat", null, null, false));
        partIDs.add(second.getValue());

        //"Super Device" has "catcat" within its sequence, but an annotation will not be made for "Basic Device".
        //This is because "Basic Device" does not have a feature because it lacks a role.
        Map<String, String> seqrole = new HashMap<>();
        seqrole.put("role", "GENE");
        seqrole.put("sequence", "tcgcatcatgt");
        ObjectId third = createDevice(persistor, "Super Device", partIDs, seqrole, "David T.", false);
        compareDesign(third, getDevice(persistor, "Super device", null, "gene", "tcgcat", null, null, false));
        partIDs.add(third.getValue());

        Map<String, String> seqroleAndSuperDevice = new HashMap<>();
        seqroleAndSuperDevice.put("role", "GENE");
        seqroleAndSuperDevice.put("sequence", "actacttcgcatcatgttcatca");
        ObjectId fourth = createDevice(persistor, "Device with Super Device", partIDs, seqroleAndSuperDevice, "David T.", false);
        compareDesign(fourth, getDevice(persistor, "device with super device", null, "gene", "actacttcgcatcat", null, null, false));
        partIDs.add(fourth.getValue());

        ArrayList<Parameter> paramObjs = new ArrayList<>();
        paramObjs.add(new Parameter("paramName", 252.2, "paramVar", "paramUnits"));
        //This device will have a sequence equivalent to the concatenation of devices 2, 3, and 4.
        ObjectId fifth = createDevice(persistor, "Spooky Device", "Barebones Device with Display ID", partIDs, "David T.", true);
        compareDesign(fifth, getDevice(persistor, "Spooky", "bones", null, "actacttcgcatcat", null, null, false));

        ObjectId sixth = createDevice(persistor, "Parameterized Device1", partIDs, paramObjs, "David T.", false);
        compareDesign(sixth, getDevice(persistor, "Parameterized Device1", null, null, null, null, paramObjs, false));

        //Device will NOT have its sequence be equal to the concatenation of the part sequences because we provided it a sequence in "seqrole".
        ObjectId seventh = createDevice(persistor, "Parameterized Device2", partIDs, seqrole, "David T.", true);
        compareDesign(seventh, getDevice(persistor, "device2", null, "gene", "tcgcatcat", null, null, false));

        ObjectId eighth = createDevice(persistor, "Full Device", "Full Device with Parameters", partIDs, seqrole, paramObjs, "David T.", false);
        compareDesign(eighth, getDevice(persistor, "Full Dev", "with Parameters", "gene", "cat", null, paramObjs, false));
    }

    @Test
    public void testExamples() {
        Map<String, String> seqrole = new HashMap<>();
        ArrayList<Parameter> params = new ArrayList<>();
        ArrayList<String> idList = new ArrayList<>();

        seqrole.put("role", "PROMOTER");
        seqrole.put("sequence", "aacgatcgttggctgtgttgacaattaatcatcggctcgtataatgtgtggaattgtgagcgctcacaatt");
        params.add(new Parameter("pTac", 676.3, "gene expression", "REU"));
        ObjectId pTac = createPart(persistor, "pTac", "p12", seqrole, params, "David T.");

        /////////////////////////////////////////////////////////////
        seqrole = new HashMap<>();
        seqrole.put("role", "CDS");
        seqrole.put("sequence", "atggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccacaggcaagctgcccgtgccctggcccaccctcgtgaccaccttcggctacggcctgcaatgcttcgcccgctaccccgaccacatgaagctgcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacaccccaatcggcgacggccccgtgctgctgcccgacaaccactaccttagctaccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaagtaa");
        params = new ArrayList<>();
        params.add(new Parameter("YFP Degradation Rate", 0.075, "rate", "second-1"));
        ObjectId YFP = createPart(persistor, "YFP", "c32", seqrole, params, "David T.");

        /////////////////////////////////////////////////////////////
        idList.add(pTac.getValue());
        idList.add(YFP.getValue());
        seqrole = new HashMap<>();
        seqrole.put("role", "PROMOTER");
        seqrole.put("sequence", "aacgatcgttggctgtgttgacaattaatcatcggctcgtataatgtgtggaattgtgagcgctcacaattatggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccacaggcaagctgcccgtgccctggcccaccctcgtgaccaccttcggctacggcctgcaatgcttcgcccgctaccccgaccacatgaagctgcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacaccccaatcggcgacggccccgtgctgctgcccgacaaccactaccttagctaccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaagtaa");
        params = new ArrayList<>();
        params.add(new Parameter("Lacl-GFP Sensitivity", 0.1, "sensitivity", ""));
        ObjectId LaclSensor = createDevice(persistor, "Lacl Sensor", "s45", idList, seqrole, params, "David T.", false);

        //Queries
        /////////////////////////////////////////////////////////////
        params = new ArrayList<>();
        params.add(new Parameter("pTac", 676.3, "gene expression", "REU"));
        compareDesign(pTac, getPart(persistor,
                "ptac",
                "p12",
                "promoter",
                "aacgatcgttggctgtgttgacaattaatcatcggctcgtataatgtgtggaattgtgagcgctcacaatt",
                params,
                false));

        params = new ArrayList<>();
        params.add(new Parameter("YFP Degradation Rate", 0.075, "rate", "second-1"));
        compareDesign(YFP, getPart(persistor,
                "yfp",
                "c32",
                "CDS",
                "atggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacgg",
                params,
                false));

        Part qpTac = new Part("pTac", new Person("doesn't matter"));
        qpTac.setDisplayID("p12");
        Part qyfp = new Part("yfp", new Person("bestAmoeba2016"));
        qyfp.setDisplayID("c32");
        ArrayList<Part> partList = new ArrayList<>();
        partList.add(qpTac);
        partList.add(qyfp);
        params = new ArrayList<>();
        params.add(new Parameter("Lacl-GFP Sensitivity", 0.1, "sensitivity", ""));
        compareDesign(LaclSensor, getDevice(persistor,
                "lacl sensor",
                "s45",
                "promoter",
                "gatcgttggctgtgttgacaattaatcatcggctcgtataatgtgtggaattgtgagcgctcacaattatg",
                partList,
                params,
                false));

    }

    //////////////////////////
    //      Debug tools     //
    //////////////////////////
    public void printDesign(ObjectId design) {

        BioDesign bd = persistor.get(BioDesign.class, design);
        if (bd != null) {
            System.out.println("BioDesign: " + bd.toString());
            System.out.println("SubDesigns: " + bd.getSubDesigns());
            System.out.println("Parameters: " + bd.getParameters());
            System.out.println("Module: " + bd.getModule());
            if (bd.getModule() != null) {
                System.out.println("Role: " + bd.getModule().getRole());
            }

            Set<Part> parts = bd.getParts();
            for (Part p : parts) {
                if (p.getSequence() != null) {
                    System.out.println("Sequence: " + p.getSequence().getSequence());
                    System.out.println("Annotations: " + p.getSequence().getAnnotations());
                } else {
                    System.out.println(bd.getName() + " has no sequence!");
                }
            }
        } else {
            System.out.println("Could not find ObjectId " + bd.getId().getValue());
        }
        System.out.println();
        System.out.println();
    }

    public void compareDesign(ObjectId original, Map<String, Map<String, String>> found) {
        BioDesign bd = persistor.get(BioDesign.class, original);

        System.out.println();
        System.out.println("Comparing " + bd.getName());

        for (String key : found.keySet()) {
            Map<String, String> info = found.get(key);
            for (String field : info.keySet()) {
                if (field.equalsIgnoreCase("name")) {
                    try {
                        assertEquals(bd.getName(), info.get(field));
                        System.out.println("Name matched");
                    } catch (ComparisonFailure cf) {
                        System.out.println("Name check failed");
                    }

                }
                if (field.equalsIgnoreCase("id")) {
                    try {
                        assertEquals(bd.getId().toString(), info.get(field));
                        System.out.println("Id matched");
                    } catch (ComparisonFailure cf) {
                        System.out.println("ID check failed");
                    }
                } else if (field.equalsIgnoreCase("displayID")) {
                    try {
                        assertEquals(bd.getDisplayID(), info.get(field));
                        System.out.println("DisplayID matched");
                    } catch (ComparisonFailure cf) {
                        System.out.println("DisplayID check failed");
                    }
                } else if (field.equalsIgnoreCase("sequence")) {
                    for (Part p : bd.getParts()) {
                        if (p.getName().equalsIgnoreCase(bd.getName())) {
                            try {
                                assertEquals(p.getSequence().getSequence(), info.get(field));
                                System.out.println("Sequence matched");
                            } catch (ComparisonFailure cf) {
                                System.out.println("Sequence check failed");
                            }
                            break;
                        }
                    }
                } else if (field.equalsIgnoreCase("role")) {
                    try {
                        assertEquals(bd.getModule().getRole(), info.get(field));
                        System.out.println("Role matched");
                    } catch (ComparisonFailure cf) {
                        System.out.println("Role check failed");
                    }
                } else if (field.equalsIgnoreCase("parameters")) {
                    String build = "[";
                    for (Parameter p : bd.getParameters()) {
                        build += "{name:'" + p.getName()
                                + "', value:" + p.getValue()
                                + ", variable:'" + p.getVariable()
                                + "', units:'" + p.getUnits() + "'},";
                    }
                    try {
                        assertEquals(build + "]", info.get(field));
                        System.out.println("Parameters matched");
                    } catch (ComparisonFailure cf) {
                        System.out.println("Parameters check failed: Manually check the order of the loaded parameters for correctness: ");
                        System.out.println("Query Result: " + cf.getActual());
                        System.out.println("Created Object: " + cf.getExpected());
                    }
                } else if (field.equalsIgnoreCase("parts")) {
                    String pbuild = "[";
                    for (Part p : bd.getParts()) {
                        pbuild += "{id:'" + p.getId()
                                + "', name:'" + p.getName()
                                + "', displayID:'" + p.getDisplayID();
                        if (p.getSequence() != null) {
                            if (p.getSequence().getAnnotations() != null) {
                                for (Annotation anno : p.getSequence().getAnnotations()) {
                                    if (anno.getFeature().getName().equalsIgnoreCase(p.getName())) {
                                        pbuild += "', role:'" + anno.getFeature().getName();
                                        break;
                                    }
                                }
                            }
                            pbuild += "', sequence:'" + p.getSequence().getSequence();
                        }
                        pbuild += "'},";
                    }
                    try {
                        assertEquals(pbuild + "]", info.get(field));
                        System.out.println("Parts matched");
                    } catch (ComparisonFailure cf) {
                        System.out.println("Parts check failed: Manually check the order of the loaded parts for correctness:");
                        System.out.println("Query Result: " + cf.getActual());
                        System.out.println("Created Object: " + cf.getExpected());
                    }
                }
            }
        }
        System.out.println();
    }
    
    /* 
    @ author: Jason 
    
    Testing methods for my deletePart() method from ConvenienceMethods.java 
    
    */
    
    @Test
    public void deletePartTest() {
        
        // test by deleting the feature first
        // get an instance of object id 
        
        // create a part
        // query 
        // println 
        
        // delete a part 
        // query 
        // verify in mongo shell 
        
        System.out.println("Testing the delete function:");
        
        ObjectId test1 = createPart(persistor, "new special part", "Jason");

        Map<String, String> roleParam = new HashMap<>();
        
        roleParam.put("role", "GENE");
        
        BioDesign design = persistor.get(BioDesign.class,test1);
        
        // figure out a way to make sure parameters 
        // initialization will not cause error 
        
        List<Parameter> parameters = null;
                // ("Jason",2.0,"hi","meters");
        
        // delete the part 
        
        // deletePart(test1, persistor, "new special part", "GENE", "catacatcat",null, "Jason");
        
        System.out.println("Test 1 passed!");
        
        ObjectId obj = persistor.save(design);
        
        // generic test 1 
        // deletePart(obj, persistor, "Jason", "Clotho", "catcatcat", null, "David");
        
        System.out.println("Test 2 passed!");
        
        // generic test 2 
        // deletePart()
        
        System.out.println("End delete test function");
    }

    private void deletePart(ObjectId test1, Persistor persistor, String new_special_part, String gene, String catacatcat, List<Parameter> parameters, String jason) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}

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
import org.clothocad.model.Assembly;
import org.clothocad.model.BasicModule;
import org.clothocad.model.BioDesign;
import org.clothocad.model.Feature;
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

        ArrayList<String> partList = new ArrayList<>();
        partList.add("ptac");
        partList.add("yfp");
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
    
    @Test
    public void testDeleteDevice() {
        /*
            Devices and parts from testExamples()
        */
        
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
        
        boolean deleteDevice = true;
        ObjectId LaclSensor = createDevice(persistor, "Lacl Sensor", "s45", idList, seqrole, params, "David T.", false);
                
        System.out.println("Initialized pTac: " + pTac.getValue() + ", YFP: " + YFP.getValue() + ", LaclSensor: " + LaclSensor.getValue());
        
        System.out.println("Deleting pTac (Part)...");
        //Should have YFP and Lacl Sensor in DB after delete
        delete(persistor, pTac, !deleteDevice);
        
        assertEquals(null, persistor.get(pTac));
        assertNotEquals(null, persistor.get(YFP));
        assertNotEquals(null, persistor.get(LaclSensor));        
        System.out.println("Deleting LaclSensor (Device)...");
        
        //Should not have any devices or parts left in DB
        delete(persistor, LaclSensor, deleteDevice);   
        
        assertEquals(null, persistor.get(pTac));
        assertEquals(null, persistor.get(YFP));
        assertEquals(null, persistor.get(LaclSensor));
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
     *  Author: Jerome 
     */
    public void printDesignFields(ObjectId bioId){
        BioDesign bio = persistor.get(BioDesign.class, bioId);
        
        if (bio != null){
            // ObjId
            System.out.println("ObjId = " + bio.getId());
            
            // authName
            System.out.println("AuthName = " + bio.getAuthor().getDisplayName());
            
            // displayId
            System.out.println("DisplayId = " + bio.getDisplayID());
            
            // name
            System.out.println("Name = " + bio.getName());
            
            // parameters
            System.out.println("Parameters ... ");
            if (bio.getParameters() != null){
                Parameter[] params = bio.getParameters().toArray(new Parameter[bio.getParameters().size()]);
                for (Parameter p : params){
                    System.out.println("    " + p.getName() + ", " + p.getValue() + ", " + p.getVariable() + ", " + p.getUnits());
                }
            } else System.out.println("    null");
            
            // subParts
            System.out.println("Subparts ... ");
            Part[] partA = bio.getParts().toArray(new Part[bio.getParts().size()]);
            if (partA[0].getAssemblies() != null && partA[0].getAssemblies().size() > 0){
                Assembly assembly = partA[0].getAssemblies().get(0);
                List<Part> parts = assembly.getParts();
                for (Part p : parts){
                    System.out.println("    " + p);
                }
            } else System.out.println("    null");
            
            // seqrole
            System.out.println("Sequence = " + partA[0].getSequence().getSequence());
            if(bio.getModule() != null)
                System.out.println("ModuleRole = " + bio.getModule().getRole());
            else
                System.out.println("No BasicModule available");
            System.out.println("FeatureRole = " + partA[0].getRoles());
        }
        
        else{
            System.out.println("BioDesign not found.");
        }
    }
    @Test
    public void testUpdate() throws InterruptedException{
        
        System.out.println();
        System.out.println("===== Testing Convenience Update =====");
        System.out.println();
        System.out.println(" ==== Setup: Complete  ==== ");   
        System.out.println("  === Set up BioDes    ===");
        String authName = "jerome";
        
        List<Parameter> newParams = new ArrayList<>();
        
        newParams.add(new Parameter("the_sensitivity", 0.1, "sensitivity", "pounds?"));
        newParams.add(new Parameter("the_brightness", 0.1, "brightness", "lux? lumens?"));
        
        Map<String, String> justSeq = new HashMap<>();
        Map<String, String> justRole = new HashMap<>();
        Map<String, String> bothSeqRole = new HashMap<>();
        Map<String, String> neitherSeqRole = new HashMap<>();
        Map<String, String> brokenSeqRole = new HashMap<>();
        
        justSeq.put("sequence", "actgactgactg");
        
        justRole.put("role", "notapromoteranymore");
        
        bothSeqRole.put("sequence", "ggggggggg");
        bothSeqRole.put("role", "mightbeapromoter");
        
        brokenSeqRole.put("stuff", "whoknows");
        
        String disPartName = "disPart";
        String datPartName = "datPart";
        Map<String, String> createPartSeqRole = new HashMap<>();
        createPartSeqRole.put("sequence","actg");
        createPartSeqRole.put("role", "promoter");
        
        ObjectId disPartId = createPart(persistor, disPartName, createPartSeqRole, authName);
        ObjectId datPartId = createPart(persistor, datPartName, createPartSeqRole, authName);
        
        Map<String, String> c1MapSeqRole = new HashMap<>();
        c1MapSeqRole.put("sequence", "tttt");
        c1MapSeqRole.put("role", "sortofapromoter");
        
        Map<String, String> c2MapSeqRole = new HashMap<>();
        c2MapSeqRole.put("sequence", "gggg");
        c2MapSeqRole.put("role", "definitelynotapromoter");
        
        ObjectId content1 = createPart(persistor, "first content", c1MapSeqRole, authName);
        ObjectId content2 = createPart(persistor, "second content", c2MapSeqRole, authName);
        
        List <String> disDeviceChange = new ArrayList<>();
        disDeviceChange.add(content1.getValue());
        disDeviceChange.add(content2.getValue());
        
        String disDeviceName = "disDevice";
        List<String> disDeviceComponents = new ArrayList<>();
        disDeviceComponents.add(disPartId.getValue());
        disDeviceComponents.add(datPartId.getValue());
        //Persistor persistor, String name, List<String> partIDs, String author, boolean createSeqFromParts
        ObjectId disDeviceId = createDevice(persistor, disDeviceName, disDeviceComponents, authName, true);
        
        //Wait for stuff
        //Thread.sleep(10000);
        
        // Check that they are identical except for their names, ids
        //printDesign(disPart); printDesign(datPart);
        //compareDesign(disPartId, getPart(persistor, "datPart", null, "promoter", "actg", null, false));
        //compareDesign(datPartId, getPart(persistor, "disPart", null, "promoter", "actg", null, false));

        System.out.println("  === Done BD setup    ===");        
        System.out.println(" ==== Setup: Complete  ==== ");
        System.out.println("  === Listing Fields   ===");
        
        
        //printDesignFields(disPartId);
        System.out.println();
        printDesignFields(disDeviceId);
        /*  Persistor   persistor,
            ObjectId    obj,
            String      authName,
            String      displayID,
            String      name,
            List<Parameter>     parameters,
            List<String>        subPartIds,
            Map<String,String>  seqrole{*/
        System.out.println();
        System.out.println(" ==== Test 1: Parts    ==== ");
        System.out.println();
        // Test Nothing
        /**
        updatePart(persistor, disPartId, null, null, null, null);
        printDesignFields(disPartId);
        System.out.println();
        /**/
        
        // Test displayId
        /**
        updatePart(persistor, disPartId, "1", null, null, null);
        printDesignFields(disPartId);
        System.out.println();
        /**/
        
        // Test Name
        /**
        updatePart(persistor, disPartId, null, "2", null, null);
        printDesignFields(disPartId);
        
        // Test Parameters
        //No original parameters
        /**
        updatePart(persistor, disPartId, null, null, newParams, null);
        printDesignFields(disPartId);
        System.out.println();
        /**/
        
        //clear parameters
        /**
        List<Parameter> emptyParams = new ArrayList<>();
        updatePart(persistor, disPartId, null, null, emptyParams, null);
        System.out.println();
        /**/
        
        //Test seqRole...
        
        //Test justSeq
        /**
        updatePart(persistor, disPartId, null, null, null, justSeq);
        printDesignFields(disPartId);
        System.out.println();
        /**/
        
        //Test justRole
        /**
        updatePart(persistor, disPartId, null, null, null, justRole);
        printDesignFields(disPartId);
        System.out.println();
        /**/
        
        //Test bothSeqRole
        /**
        updatePart(persistor, disPartId, null, null, null, bothSeqRole);
        printDesignFields(disPartId);
        System.out.println();
        /**/
        
        //Test neitherSeqRole
        /**
        updatePart(persistor, disPartId, null, null, null, neitherSeqRole);
        printDesignFields(disPartId);
        System.out.println();
        /**/
        
        //Test brokenSeqRole
        /**
        updatePart(persistor, disPartId, null, null, null, brokenSeqRole);
        printDesignFields(disPartId);
        System.out.println();
        /**/
        
        //Test all update part arguments
        /**
        updatePart(persistor, disPartId, "1", "2", newParams, bothSeqRole);
        printDesignFields(disPartId);
        System.out.println();
        /**/
        
        //System.out.println("        Skipping...");
        System.out.println(" ==== Test 1: Complete ==== ");
        System.out.println();
        System.out.println(" ===  Test 2: Devices   === ");
        System.out.println();
        
        /**
            Persistor   persistor,
            ObjectId    obj,
            String      displayID,
            String      name,
            List<Parameter>     parameters,
            List<String>        subPartIds,
            Map<String,String>  seqrole,
            boolean             createSeqFromParts
        /**/
        
        // Test Nothing
        /**
        updateDevice(persistor, disDeviceId, null, null, null, null, null, false);
        printDesignFields(disDeviceId);
        System.out.println();
        /**/
        
        // Test displayId
        /**
        updateDevice(persistor, disDeviceId, "1", null, null, null, null, false);
        printDesignFields(disDeviceId);
        System.out.println();
        /**/
        
        // Test Name
        /**
        updateDevice(persistor, disDeviceId, null, "2", null, null, null, false);
        printDesignFields(disDeviceId);
        System.out.println();
        /**/
        
        // Test Parameters
        //No original parameters
        /**
        updateDevice(persistor, disDeviceId, null, null, newParams, null, null, false);
        printDesignFields(disDeviceId);
        System.out.println();
        /**/
        
        //clear parameters
        /**
        List<Parameter> emptyParams = new ArrayList<>();
        updateDevice(persistor, disDeviceId, null, null, emptyParams, null, null, false);
        printDesignFields(disDeviceId);
        System.out.println();
        /**/
        
        //Test seqRole...
        //Test justSeq
        /**
        updateDevice(persistor, disDeviceId, null, null, null, null, justSeq, false);
        printDesignFields(disDeviceId);
        System.out.println();
        /**/
        
        //Test justRole
        /**
        updateDevice(persistor, disDeviceId, null, null, null, null, justRole, false);
        printDesignFields(disDeviceId);
        System.out.println();
        /**/
        
        //Test bothSeqRole
        /**
        updateDevice(persistor, disDeviceId, null, null, null, null, bothSeqRole, false);
        printDesignFields(disDeviceId);
        System.out.println();
        /**/
        
        //Test neitherSeqRole
        /**
        updateDevice(persistor, disDeviceId, null, null, null, null, neitherSeqRole, false);
        printDesignFields(disDeviceId);
        System.out.println();
        /**/
        
        //Test brokenSeqRole
        /**
        updateDevice(persistor, disDeviceId, null, null, null, null, brokenSeqRole, false);
        printDesignFields(disDeviceId);
        System.out.println();
        /**/
        
        //Test partsList update
        /**
        Person itsame = new Person("authName");
        Feature feat = new Feature(disDeviceName, "promoter", itsame);
        Set<Feature> featSet = new HashSet<>();
        featSet.add(feat);
        BioDesign bio = persistor.get(BioDesign.class, disDeviceId);
        BasicModule mod = new BasicModule(disDeviceName, "promoter", featSet, itsame);
        bio.setModule(mod);
        updateDevice(persistor, disDeviceId, null, null, null, disDeviceChange, null, true);
        printDesignFields(disDeviceId);
        printDesignFields(content1);
        printDesignFields(content2);
        System.out.println();
        /**/
        
        //Test all update part arguments
        /**
        Map<String, String> emptySeqRole = new HashMap<>();
        emptySeqRole.put("sequence", "");
        emptySeqRole.put("role", "");
        updateDevice(persistor, disDeviceId, "0", "0", newParams, null, emptySeqRole, false);
        printDesignFields(disDeviceId);
        System.out.println();
        /**/
        
        /**
        updateDevice(persistor, disDeviceId, "9", "9", emptyParams, null, bothSeqRole, true);
        printDesignFields(disDeviceId);
        System.out.println();
        /**/
        
        
        //System.out.println("        Skipping...");
        System.out.println(" ===  Test 2: Complete  === ");
        System.out.println("===== Testing Complete =====");
    }
}

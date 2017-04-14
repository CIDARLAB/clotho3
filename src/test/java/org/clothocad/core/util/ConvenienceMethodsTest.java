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
import java.util.Map;
import java.util.Set;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.model.BioDesign;
import org.clothocad.model.Parameter;
import org.clothocad.model.Part;
import org.clothocad.webserver.jetty.ConvenienceMethods;
import static org.clothocad.webserver.jetty.ConvenienceMethods.createDevice;
import static org.clothocad.webserver.jetty.ConvenienceMethods.createPart;
import org.junit.Test;

/**
 *
 * @author David
 */
public class ConvenienceMethodsTest extends AuthorizedShiroTest {

    private static Persistor persistor;

    public ConvenienceMethodsTest() {
        persistor = injector.getInstance(Persistor.class);
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
        paramObj.add(new Parameter("pTac", 626.3,"gene expression", "REU"));
        ObjectId fifth = createPart(persistor, "PartWithName", "Hello I am Name", "David");
        
        ObjectId sixth = createPart(persistor, "AnotherPart", "Another worthy part", paramObj, "David");
        
        ObjectId seventh = createPart(persistor, "FullPart", "Full Part with all Parameters", seqrole, paramObj, "David");
        
        testDesign(first);
        testDesign(second);
        testDesign(third);
        testDesign(fourth);
        testDesign(fifth);
        testDesign(sixth);
        testDesign(seventh);
    }

    @Test
    public void testCreateDevice() {

        ArrayList<String> partIDs = new ArrayList<>();

        ObjectId first = createDevice(persistor, "Barebones Device", partIDs, "David T.", false);

        partIDs.add(first.getValue());

        Map<String, String> sequence = new HashMap<>();
        sequence.put("sequence", "catcat");
        ObjectId second = createDevice(persistor, "Basic Device", partIDs, sequence, "David T.", false);

        partIDs.add(second.getValue());

        //"Super Device" has "catcat" within its sequence, but an annotation will not be made for "Basic Device".
        //This is because "Basic Device" does not have a feature because it lacks a role.
        Map<String, String> seqrole = new HashMap<>();
        seqrole.put("role", "GENE");
        seqrole.put("sequence", "tcgcatcatgt");
        ObjectId third = createDevice(persistor, "Super Device", partIDs, seqrole, "David T.", false);

        partIDs.add(third.getValue());

        Map<String, String> seqroleAndSuperDevice = new HashMap<>();
        seqroleAndSuperDevice.put("role", "GENE");
        seqroleAndSuperDevice.put("sequence", "actacttcgcatcatgttcatca");
        ObjectId fourth = createDevice(persistor, "Device with Super Device", partIDs, seqroleAndSuperDevice, "David T.", false);

        partIDs.add(fourth.getValue());

        ArrayList<Parameter> paramObjs = new ArrayList<>();
        paramObjs.add(new Parameter("paramName", 252.2, "paramVar", "paramUnits"));

        //This device will have a sequence equivalent to the concatenation of devices 2, 3, and 4.
        ObjectId fifth = createDevice(persistor, "Spooky Device", partIDs, "Barebones Device with Display ID", "David T.", true);

        ObjectId sixth = createDevice(persistor, "Parameterized Device1", partIDs, paramObjs, "David T.", false);

        //Device will NOT have its sequence be equal to the concatenation of the part sequences because we provided it a sequence in "seqrole".
        ObjectId seventh = createDevice(persistor, "Parameterized Device2", partIDs, seqrole, "David T.", true);

        ObjectId eighth = createDevice(persistor, "Full Device", partIDs, seqrole, paramObjs, "Full Device with Parameters", "David T.", false);

        testDesign(first);
        testDesign(second);
        testDesign(third);
        testDesign(fourth);
        testDesign(fifth);
        testDesign(sixth);
        testDesign(seventh);
        testDesign(eighth);
    }

    public void testDesign(ObjectId design) {

        BioDesign bd = persistor.get(BioDesign.class, design);
        if (bd != null) {
            System.out.println("BioDesign: " + bd.toString());
            System.out.println("SubDesigns: " + bd.getSubDesigns());
            System.out.println("Parameters: " + bd.getParameters());
            System.out.println("Module: " + bd.getModule());
            if(bd.getModule() != null)
            {
                System.out.println("Role: " + bd.getModule().getRole());
            }

            Set<Part> parts = bd.getParts();
            for (Part p : parts) {
                if (p.getSequence() != null) {
                    System.out.println("Sequence: " + p.getSequence().getSequence());
                    System.out.println("Annotations: " + p.getSequence().getAnnotations());
                }
                else
                {
                    System.out.println(bd.getName() + " has no sequence!");
                }
            }
        } else {
            System.out.println("Could not find ObjectId " + bd.getId().getValue());
        }
        System.out.println();
        System.out.println();
    }

}

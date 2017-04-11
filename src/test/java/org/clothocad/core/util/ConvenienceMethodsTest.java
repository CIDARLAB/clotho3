///*
// * To change this license header, choose License Headers in Project Properties.
// * To change this template file, choose Tools | Templates
// * and open the template in the editor.
// */
//package org.clothocad.core.util;
//
//import java.util.ArrayList;
//import java.util.HashMap;
//import java.util.HashSet;
//import java.util.Map;
//import java.util.Set;
//import org.clothocad.core.datums.ObjBase;
//import org.clothocad.core.datums.ObjectId;
//import org.clothocad.core.persistence.Persistor;
//import org.clothocad.webserver.jetty.ConvenienceMethods;
//import static org.clothocad.webserver.jetty.ConvenienceMethods.createDevice;
//import static org.clothocad.webserver.jetty.ConvenienceMethods.createPart;
//import org.junit.Test;
//
///**
// *
// * @author David
// */
//public class ConvenienceMethodsTest extends AuthorizedShiroTest {
//
//    private static Persistor persistor;
//
//    public ConvenienceMethodsTest() {
//        persistor = injector.getInstance(Persistor.class);
//    }
//
//    /*
//    BioDesign
//    name
//    
//    Part
//    name
//    
//    Sequence - Only if sequence provided
//    name
//    sequence
//    
//    Annotation - Only if sequence or role provided
//    name
//    start
//    end
//    
//    Feature - Only if role provided
//    name
//    role
//    
//    BasicModule - Only if role provided
//    name
//    role
//     */
//    
//    /*
//    Had to change function signature to be able to capture all of the optional combinations of String parameters
//    
//    __optionals__:
//    role
//    sequence
//    */
//    
//    @Test
//    public void testCreatePart(){
//        createPart(persistor,"mySpecialPart", "David");
//        
//        Map<String, String> roleParam = new HashMap<>();
//        roleParam.put("role", "PROMOTER");
//        createPart(persistor,"roleOnlyPart", roleParam, "David");
//        
//        Map<String, String> seqParam = new HashMap<>();
//        seqParam.put("sequence", "catcatcatcatcatcatcatcatcat");
//        createPart(persistor,"FunCatPart", seqParam, "David");
//        
//        Map<String, String> bothParams = new HashMap<>();
//        bothParams.put("role", "GENE");
//        bothParams.put("sequence", "tccctatcagtgatagagattgacatccctatcagtgatagagatactgagcac");
//        createPart(persistor, "R0040 Sequence",bothParams, "June Rhee");
//    }
//    
//    @Test
//    public void testCreateDevice(){
//        
//        ArrayList<String> partIDs = new ArrayList<>();
//        
//        ObjectId zero = createDevice(persistor, "DevicePls", partIDs, "David T.");
//        
//        partIDs.add(zero.getValue());
//        
//        Map<String, String> seqParam = new HashMap<>();
//        seqParam.put("sequence", "catcat");
//        ObjectId first = createDevice(persistor, "Basic Device", partIDs, seqParam, "David T.");
//        
//        partIDs.add(first.getValue());
//        
//        //"Super Device" has "catcat" within its sequence, but an annotation will not be made for "Basic Device".
//        //This is because "Basic Device" does not have a feature because it lacks a role.
//        Map<String, String> bothParams = new HashMap<>();
//        bothParams.put("role", "GENE");
//        bothParams.put("sequence", "tcgcatcatgt");
//        ObjectId second = createDevice(persistor, "Super Device", partIDs, bothParams, "David T.");
//        
//        partIDs.add(second.getValue());
//        
//        Map<String, String> bothParamsAndSuperDevice = new HashMap<>();
//        bothParamsAndSuperDevice.put("role", "GENE");
//        bothParamsAndSuperDevice.put("sequence", "actacttcgcatcatgttcatca");
//        ObjectId third = createDevice(persistor, "Device with Super Device", partIDs, bothParamsAndSuperDevice, "David T.");
//        
//        partIDs.add(third.getValue());
//        
//    }
//    
//    
//    @Test
//    public void testCreations(){
//        
//        Map<String, Object> query = new HashMap<>();
//        query.put("author", "david");
//        Iterable<ObjBase> list = persistor.findRegex(query);
//        
//        for(ObjBase each : list)
//        {
//            System.out.println(each.toString());
//        }
//    }
//    
//    
//    
//
//}

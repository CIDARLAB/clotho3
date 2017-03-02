/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.webserver.jetty.ConvenienceMethods;
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
    public void testCreations(){
        ConvenienceMethods.
        createPart(persistor,"mySpecialPart", "David");
        
        Map<String, String> roleParam = new HashMap<>();
        roleParam.put("role", "PROMOTER");
        createPart(persistor,"roleOnlyPart", roleParam, "David");
        
        Map<String, String> seqParam = new HashMap<>();
        seqParam.put("sequence", "catcatcatcatcatcatcatcatcat");
        createPart(persistor,"FunCatPart", seqParam, "David");
        
        Map<String, String> bothParams = new HashMap<>();
        bothParams.put("role", "GENE");
        bothParams.put("sequence", "tccctatcagtgatagagattgacatccctatcagtgatagagatactgagcac");
        createPart(persistor, "R0040 Sequence",bothParams, "June Rhee");
        
        Map<String, Object> query = new HashMap<>();
        query.put("author", "david");
        Iterable<ObjBase> list = persistor.findRegex(query);
        
        for(ObjBase each : list)
        {
            System.out.println(each.toString());
        }
        
    }
    
    
    

}

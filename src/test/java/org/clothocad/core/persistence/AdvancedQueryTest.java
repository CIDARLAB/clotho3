/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.clothocad.core.communication.Router;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.util.AuthorizedShiroTest;
import static org.clothocad.core.util.LumazineSynthaseExample.createDesignK542008;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 *
 * @author David
 */
public class AdvancedQueryTest extends AuthorizedShiroTest {

    /////////////////////////
    //
    //When testing this file, keep in mind that this file will interact with the db defined in ClothoTestModule.java.
    //The db is most likely 'testClotho'.
    //
    ///////////////////////////

    
    private static Persistor persistor;
    private static Router router;
    private boolean first = true;

    public AdvancedQueryTest() {
        if (first) {

            persistor = injector.getInstance(Persistor.class);
            router = injector.getInstance(Router.class);

            //this is expensive. Let's just initialize the data once please.
            createDesignK542008(persistor);
            first = false;
        }

        System.out.println();
        System.out.println();
    }

    
    //I'm new to tests, forgive me for making a single test and calling my test functions from it
    @Test
    public void runTests(){
        listAllTest();
        regexTest();
    }
    
    public void listAllTest() {

        List<Map<String, Object>> all = new ArrayList<>();
        Collection<ObjBase> list = persistor.listAll();

        for (ObjBase each : list) {
            System.out.println("ALL LIST : " + each);
        }
    }
    
    public void regexTest() {
        //
        // Take a look at a random line in createDesignK542008(), and try to apply regex to select some of the objects created
        // protip: go to this test class's constructor, click on the function, and hit Ctrl-B (or Mac equivalent) to fast travel
        //

        Map<String, Object> query = new HashMap<>();
        query.put("name", "bba.*K");
        Iterable<ObjBase> list = persistor.findRegex(query);
        for (ObjBase each : list) {
            System.out.println("REGEX NAME : " + each);
        }

        query = new HashMap<>();
        query.put("role", "promoter");
        list = persistor.findRegex(query);
        for (ObjBase each : list) {
            System.out.println("REGEX ROLE : " + each);
        }
        
        query = new HashMap<>();
        query.put("sequence", "cat.*cat.*tgggg");
        list = persistor.findRegex(query);
        for (ObjBase each : list) {
            System.out.println("REGEX SEQUENCE : " + each);
        }
    }
}

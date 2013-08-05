package org.clothocad.core.testers;

import com.google.inject.Injector;
import org.clothocad.core.layers.communication.ScriptAPI;
import org.clothocad.core.layers.communication.connection.ClientConnection;
import org.junit.Before;
import org.junit.Test;
import org.clothocad.core.layers.communication.mind.Mind;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.testers.communication.TestConnection;
import org.clothocad.core.util.FileUtils;
import org.clothocad.core.utils.TestUtils;
import org.junit.After;
import static org.junit.Assert.*;

public class MindTest {
    Mind mind = new Mind();
    
    @Before
    public void setUp() {
        System.out.println("setup");
        //setup test environment
        
        //add specific test data
    }
    
    @After
    public void tearDown() {
        System.out.println("tear");
    }

/*    @Test
    public void testJSScript() throws Exception{
        String script = FileUtils.readFile("mind_test.js");
        assertTrue(mind.runCommand(script, new ScriptAPI(mind, null, "")));
    }*/
    
    @Test
    public void testScript() {
        TestConnection connection = new TestConnection("testScriptGet");
        mind.setConnection(connection);
        Injector injector = TestUtils.getDefaultTestInjector();
        Persistor persistor = injector.getInstance(Persistor.class);
        persistor.deleteAll();
        String script = "var newobjid = clotho.create( {\"name\":\"UCM\",\"state\":\"MA\",\"schema\":\"Institution\",\"country\":\"United States of America\",\"city\":\"Baltizam\"} );\n" +
                        "var result = clotho.get(newobjid);\n" +
                        "if(result.name != \"UCM\") {\n" +
                        "    throw \"wrong name: \" + result.name;\n" +
                        "}";
        assertTrue(mind.runCommand(script, new ScriptAPI(mind, persistor, "")));
        script = "result = clotho.get('UCM');\n" +
                 "if(result.name != \"UCM\") {\n" +
                 "    throw \"wrong name: \" + result.name;\n" +
                 "}";
        assertTrue(mind.runCommand(script, new ScriptAPI(mind, persistor, "")));
        script = "wrapper = [];\n" +
                 "wrapper[0] = newobjid;\n" +
                 "result = clotho.get(wrapper);\n" +
                 "if(result.name != \"UCM\") {\n" +
                 "    throw \"wrong name: \" + result.name;\n" +
                 "}";
        assertTrue(mind.runCommand(script, new ScriptAPI(mind, persistor, "")));
        
        script = "var listy = clotho.query({\"city\" : \"Baltizam\"});\n" +
                 "if (listy.length != 1) {\n" +
                 "   throw \"wrong number of results: \" + listy.size;\n" +
                 "}\n" +
                 "var existing = listy[0];\n" +
                 "if(existing.name != \"UCM\") {\n" +
                 "    throw \"wrong name: \" + existing.name;\n" +
                 "}";
        assertTrue(mind.runCommand(script, new ScriptAPI(mind, persistor, "")));
        
        script = "var args = {};\n" +
                 "args.id = newobjid;\n" +
                 "args.city = \"Paris\";\n" +
                 "result = clotho.set(args);\n" +
                 "if(result.city != \"Paris\") {\n" +
                 "    throw \"wrong city: \" + result.city;\n" +
                 "}";
        assertTrue(mind.runCommand(script, new ScriptAPI(mind, persistor, ""))); 
        
        script = "existing.name = \"Shamoo University\"; \n" +
                 "existing.city = \"Whaletown\"; \n" +
                 "existing.state = \"NR\"; \n" +
                 "result = clotho.set(existing);\n" +
                 "if(result.city != \"Whaletown\") {\n" +
                 "    throw \"wrong city: \" + result.city;\n" +
                 "}\n" +
                 "if(result.state != \"NR\") {\n" +
                 "    throw \"wrong state: \" + result.state;\n" +
                 "}";
        assertTrue(mind.runCommand(script, new ScriptAPI(mind, persistor, ""))); 
        
        script = "var finalSet = clotho.query({\"city\" : \"Whaletown\"});\n" +
                 "if(finalSet.length!=1) {\n" +
                 "   throw \"wrong number of results: \" + finalSet.size;\n" +
                 "}";
        assertTrue(mind.runCommand(script, new ScriptAPI(mind, persistor, ""))); 
        
        script = "clotho.destroy(newobjid);\n" +
                 "finalSet = clotho.query({\"city\" : \"Whaletown\"});\n" +
                 "if(finalSet.length!=0) {\n" +
                 "   throw \"wrong number of results: \" + finalSet.size;\n" +
                 "}";
        assertTrue(mind.runCommand(script, new ScriptAPI(mind, persistor, ""))); 
    }
}

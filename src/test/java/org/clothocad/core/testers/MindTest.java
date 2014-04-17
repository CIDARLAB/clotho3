package org.clothocad.core.testers;

import com.google.inject.Injector;
import javax.script.Bindings;
import javax.script.ScriptContext;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;
import org.clothocad.core.communication.MessageOptions;
import org.clothocad.core.communication.Router;
import org.clothocad.core.execution.ScriptAPI;
import org.junit.Before;
import org.junit.Test;
import org.clothocad.core.execution.Mind;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.communication.TestConnection;
import org.clothocad.core.util.TestUtils;
import org.junit.After;

public class MindTest {
    Mind mind = new Mind();
    
    @Before
    public void setUp() {
        System.out.println("setup");
    }
    
    @After
    public void tearDown() {
        System.out.println("tear");
    }
    
    @Test
public void testScript() throws ScriptException {
        TestConnection connection = new TestConnection("testScriptGet");
        mind.setConnection(connection);
        Injector injector = TestUtils.getDefaultTestInjector();
        Persistor persistor = injector.getInstance(Persistor.class);
        Router router = injector.getInstance(Router.class);
        persistor.deleteAll();
        String script = "var newobjid = clotho.create( {\"name\":\"UCM\",\"state\":\"MA\",\"schema\":\"Institution\",\"country\":\"United States of America\",\"city\":\"Baltizam\"} );\n" +
                        "var result = clotho.get(newobjid);\n" +
                        "if(result.name != \"UCM\") {\n" +
                        "    throw \"wrong name: \" + result.name;\n" +
                        "}";
        mind.runCommand(script, new ScriptAPI(mind, persistor, router, "", new MessageOptions()));
        script = "result = clotho.get('UCM');\n" +
                 "if(result.name != \"UCM\") {\n" +
                 "    throw \"wrong name: \" + result.name;\n" +
                 "}";
        mind.runCommand(script, new ScriptAPI(mind, persistor, router, "", new MessageOptions()));
        script = "wrapper = [];\n" +
                 "wrapper[0] = newobjid;\n" +
                 "result = clotho.get(wrapper);\n" +
                 "if(result.name != \"UCM\") {\n" +
                 "    throw \"wrong name: \" + result.name;\n" +
                 "}";
        mind.runCommand(script, new ScriptAPI(mind, persistor, router, "", new MessageOptions()));
        
        script = "var listy = clotho.query({\"city\" : \"Baltizam\"});\n" +
                 "if (listy.length != 1) {\n" +
                 "   throw \"wrong number of results: \" + listy.size;\n" +
                 "}\n" +
                 "var existing = listy[0];\n" +
                 "if(existing.name != \"UCM\") {\n" +
                 "    throw \"wrong name: \" + existing.name;\n" +
                 "}";
        mind.runCommand(script, new ScriptAPI(mind, persistor, router, "", new MessageOptions()));
        
        script = "var args = {};\n" +
                 "args.id = newobjid;\n" +
                 "args.city = \"Paris\";\n" +
                 "result = clotho.set(args);\n" +
                 "if(result.city != \"Paris\") {\n" +
                 "    throw \"wrong city: \" + result.city;\n" +
                 "}";
        mind.runCommand(script, new ScriptAPI(mind, persistor, router, "", new MessageOptions())); 
        
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
        mind.runCommand(script, new ScriptAPI(mind, persistor, router, "", new MessageOptions())); 
        
        script = "var finalSet = clotho.query({\"city\" : \"Whaletown\"});\n" +
                 "if(finalSet.length!=1) {\n" +
                 "   throw \"wrong number of results: \" + finalSet.size;\n" +
                 "}";
        mind.runCommand(script, new ScriptAPI(mind, persistor, router, "", new MessageOptions())); 
        
        script = "clotho.destroy(newobjid);\n" +
                 "finalSet = clotho.query({\"city\" : \"Whaletown\"});\n" +
                 "if(finalSet.length!=0) {\n" +
                 "   throw \"wrong number of results: \" + finalSet.size;\n" +
                 "}";
        mind.runCommand(script, new ScriptAPI(mind, persistor, router, "", new MessageOptions())); 
    }
    
    @Test
    public  void adgfdasg() throws ScriptException
   {
       
       //http://www.informit.com/articles/article.aspx?p=696621&seqNum=6
       
       
       
       
       
       
       
       
       
      // Create a ScriptEngineManager that discovers all script engine
      // factories (and their associated script engines) that are visible to
      // the current thread's classloader.

      ScriptEngineManager manager = new ScriptEngineManager ();

      // Obtain a ScriptEngine that supports the JavaScript short name.

      ScriptEngine engine = manager.getEngineByName ("JavaScript");

      // Initialize the color and shape script variables.

      engine.put ("color", "red");
      engine.put ("shape", "rectangle");

      // Evaluate a script that outputs the values of these variables.

      engine.eval ("println (color); println (shape);");

      // Save the current bindings object.

      Bindings oldBindings = engine.getBindings (ScriptContext.ENGINE_SCOPE);

      // Replace the bindings with a new bindings that overrides color and
      // shape.

      Bindings newBindings = engine.createBindings ();
      newBindings.put ("color", "blue");
      engine.setBindings (newBindings, ScriptContext.ENGINE_SCOPE);
      engine.put ("shape", "triangle");

      // Evaluate the script.

      engine.eval ("println (color); println (shape);");

      // Restore the original bindings.

      engine.setBindings (oldBindings, ScriptContext.ENGINE_SCOPE);

      // Evaluate the script.

      engine.eval ("println (color); println (shape);");
      
   }
}

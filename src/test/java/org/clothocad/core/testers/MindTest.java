package org.clothocad.core.testers;

import javax.script.Bindings;
import javax.script.ScriptContext;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;
import org.junit.Before;
import org.junit.Test;
import org.clothocad.core.layers.communication.mind.Mind;
import org.clothocad.core.util.FileUtils;
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
    public void testJSScript() throws Exception{
        String script = FileUtils.readFile("mind_test.js");
        assert(mind.runCommand(script));
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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.layers.execution;

import java.util.HashMap;
import java.util.Map;
import java.util.Properties;
import javax.script.ScriptContext;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;    

/**
 *
 * @author spaige
 */
public class ClothoScriptEngineManager {
    private static ScriptEngineManager em = new ScriptEngineManager();
    
    private static Map<String, ScriptEngine> enginesByLanguage = new HashMap<String, ScriptEngine>();
    
    private static final Map<String, ScriptAPI> apis = new HashMap<String, ScriptAPI>();
    static {
        apis.put("javascript", new JavaScriptAPI());
        apis.put("python", new PythonAPI());
        //apis.put("perl", new PerlScriptAPI());
    }
    
    static {
        Properties props = System.getProperties();
        //props.setProperty("rhino.opt.level", "9");
    }

    public static ScriptEngine getEngineByName(String name){
        if (enginesByLanguage.containsKey(name)){
            return enginesByLanguage.get(name);
        }
        ScriptEngine engine = em.getEngineByName(name);
        
        injectAPI(engine, name); //inject on Global level
        enginesByLanguage.put(name, engine);
        return engine;
    }
    
    private static void injectAPI(ScriptEngine engine , String name){
        engine.getContext().getBindings(ScriptContext.GLOBAL_SCOPE).put("clotho", apis.get(name));
        apis.get(name).setEngine(engine);
    }
}

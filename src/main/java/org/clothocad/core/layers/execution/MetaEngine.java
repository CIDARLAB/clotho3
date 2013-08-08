/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.layers.execution;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.script.Invocable;
import javax.script.ScriptContext;
import javax.script.ScriptEngine;
import javax.script.ScriptException;
import javax.script.SimpleScriptContext;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.layers.communication.ScriptAPI;
import static org.clothocad.core.layers.execution.JavaScriptScript.engine;
import sun.org.mozilla.javascript.RhinoException;

/**
 *
 * @author spaige
 */
public class MetaEngine {
    
    private Map<Language,ScriptContext> contexts = new HashMap<>();
    private Map<Language,ScriptEngine> engines = new HashMap<>();
    
    public Object eval(String script, Language language, ScriptAPI api) throws ScriptException{
        ScriptContext context = getContext(language);
        ScriptEngine engine = getEngine(language);
        
        injectAPI(api, context);
        try{
            return engine.eval(script, context);          
        }
        //XXX: ew, referring to a concrete engine implementation!
        catch (RhinoException e){
            //is this terrible?
            throw new ScriptException(e);
        }
    }

    private ScriptContext getContext(Language language) {
        if (!contexts.containsKey(language)) contexts.put(language, new SimpleScriptContext());
        return contexts.get(language);
    }

    private ScriptEngine getEngine(Language language) {
        if (!engines.containsKey(language)) engines.put(language, ClothoScriptEngineManager.getEngineByLanguage(language));
        if (engines.get(language) == null) throw new EngineNotFoundException();
        return engines.get(language);
    }

    private void injectAPI(ScriptAPI api, ScriptContext context) {
        context.setAttribute("clotho", api, ScriptContext.ENGINE_SCOPE);
    }

    public Object get(String token) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    //XXX: de-js-ify this!
    public Object invoke(String code, String name, List args) throws ScriptException {
        try {
            ScriptEngine engine = getEngine(Language.JAVASCRIPT);
            engine.eval("var " + name+ " = " + code);
            Invocable invocable = (Invocable) engine;
            return invocable.invokeFunction(name, args);
        } catch (NoSuchMethodException ex) {
            throw new ScriptException(ex);
        }
    }
    
}

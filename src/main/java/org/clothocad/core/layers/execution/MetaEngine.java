/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.layers.execution;

import java.util.HashMap;
import java.util.Map;
import javax.script.ScriptContext;
import javax.script.ScriptEngine;
import javax.script.ScriptException;
import javax.script.SimpleScriptContext;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.layers.communication.ScriptAPI;

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
        
        return engine.eval(script, context);
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
    
}

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.layers.execution;

import com.github.jmkgreen.morphia.annotations.NotSaved;
import com.github.jmkgreen.morphia.annotations.PostLoad;
import com.github.jmkgreen.morphia.annotations.PreLoad;
import com.github.jmkgreen.morphia.annotations.PrePersist;
import com.mongodb.BasicDBObject;
import com.mongodb.DBObject;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.script.Bindings;
import javax.script.ScriptContext;
import javax.script.ScriptEngine;
import javax.script.ScriptException;
import javax.script.SimpleBindings;
import javax.script.SimpleScriptContext;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.layers.communication.ScriptAPI;
import org.mozilla.javascript.RhinoException;

/**
 *
 * @author spaige
 */
public class MetaEngine {
    
    @NotSaved
    private Map<Language,ScriptContext> contexts = new HashMap<>();
    
    @NotSaved
    private Map<Language,ScriptEngine> engines = new HashMap<>();
    
    private Map<Language, Map<String, Object>> bindings = new HashMap<>();
    

    
    @PreLoad
    public DBObject preLoad(DBObject obj){
        Map map = ((BasicDBObject) obj.get("bindings")).toMap();
        for (Object olang : ((BasicDBObject) obj.get("bindings")).toMap().keySet()){
            Language language = Language.valueOf(olang.toString());
            BasicDBObject b = (BasicDBObject) map.get(olang);
            bindings.put(language, b.toMap());
        }
        obj.removeField("bindings");
        return obj;
    }
    
    @PostLoad
    public void postLoad(){
        for (Language key : bindings.keySet()){
            ScriptContext context = getContext(key);
            context.setBindings(new SimpleBindings(bindings.get(key)), ScriptContext.ENGINE_SCOPE);
        }
    }
    
    @PrePersist
    public void prePersist(){
        for (Language key : contexts.keySet()){
            Bindings bindings = contexts.get(key).getBindings(ScriptContext.ENGINE_SCOPE);
            bindings.remove("clotho");
            this.bindings.put(key, bindings);
        }
    }
    
    public Object eval(String script, Language language, ScriptAPI api) throws ScriptException{
        ScriptContext context = getContext(language);
        HackEngine engine = getEngine(language);
        
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

    private HackEngine getEngine(Language language) {
        if (language == Language.JAVASCRIPT) return new JavaScriptEngine();
        
        if (!engines.containsKey(language)) engines.put(language, ClothoScriptEngineManager.getEngineByLanguage(language));
        if (engines.get(language) == null) throw new EngineNotFoundException();
        return new WrappedScriptEngine(engines.get(language));
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
            HackEngine engine = getEngine(Language.JAVASCRIPT);
            engine.eval("var " + name+ " = " + code);
            return engine.invokeFunction(name, args.toArray());
        } catch (NoSuchMethodException ex) {
            throw new ScriptException(ex);
        }
    }
    
}

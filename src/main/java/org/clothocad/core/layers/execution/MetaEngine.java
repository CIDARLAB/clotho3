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
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.script.Bindings;
import javax.script.ScriptContext;
import javax.script.ScriptEngine;
import javax.script.ScriptException;
import javax.script.SimpleBindings;
import javax.script.SimpleScriptContext;
import org.bson.types.ObjectId;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.Module;
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
        
        if (!engines.containsKey(language)) {
            engines.put(language, ClothoScriptEngineManager.getEngineByLanguage(language));
            ScriptEngine engine = engines.get(language);
            if (engine == null) throw new EngineNotFoundException();
            engine.setContext(getContext(language));
        }
        return new WrappedScriptEngine(engines.get(language));
    }

    private void injectAPI(ScriptAPI api, ScriptContext context) {
        context.setAttribute("clotho", api, ScriptContext.ENGINE_SCOPE);
    }

    public Object get(String token) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    
    public Object invoke(String functionCode, String name, List args, ScriptAPI api) throws ScriptException {
        return invoke(functionCode, "", name, args, api);
    }
    public static String idToName(ObjectId id){
        return "f"+id.toString();
    }
    
    public void loadAsGlobal(Module module, ScriptAPI api) throws ScriptException{
        String name = "f"+module.getUUID().toString();
        //TODO: check last saved date
        //if (getEngine(module.getLanguage()).getContext().getBindings(ScriptContext.ENGINE_SCOPE).containsKey(name)){
        //    return;
        //}
        HackEngine engine = getEngine(module.getLanguage());
        
        engine.eval("var " + name + " = " +module.getCodeToLoad());
    }
    
    public Object invoke(Function module, List args, ScriptAPI api) throws ScriptException, NoSuchMethodException{
       HackEngine engine = getEngine(module.getLanguage());
       injectAPI(api, engine.getContext());
       loadAsGlobal(module, api);
       
       return engine.invokeFunction(idToName(module.getUUID()), args.toArray());
    }
    
    
    public Object invoke(Module module, String methodName, List args, ScriptAPI api) throws ScriptException, NoSuchMethodException{
       HackEngine engine = getEngine(module.getLanguage());
       injectAPI(api, engine.getContext());
       loadAsGlobal(module, api);
    
       Object thiz = engine.getContext().getBindings(ScriptContext.ENGINE_SCOPE).get(idToName(module.getUUID()));
       return engine.invokeMethod(thiz, methodName, args);
    }
    
    //XXX: de-js-ify this!
    public Object invoke(String functionCode, String setupCode, String name, List args, ScriptAPI api) throws ScriptException {
        try {
            HackEngine engine = getEngine(Language.JAVASCRIPT);
            ScriptContext context = engine.getContext();
            injectAPI(api,context);
            
            //eval dependencies
            //XXX: this is terrible terrible, change to module pattern that returns a function
            engine.eval(setupCode);
            
            engine.eval("var " + name+ " = " + functionCode);
            
            
            return engine.invokeFunction(name, args.toArray());
            
        } catch (NoSuchMethodException ex) {
            throw new ScriptException(ex);
        }
    }
    
    protected String generateDependencies(Language language, Collection<ObjectId> dependencies){
        switch (language){
            case JAVASCRIPT:
                return new JavaScriptScript().generateImports(dependencies);
            default:
                throw new UnsupportedOperationException();
        }
    }
}

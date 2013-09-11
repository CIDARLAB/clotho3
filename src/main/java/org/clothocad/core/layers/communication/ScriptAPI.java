/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.layers.communication;

import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.script.ScriptException;
import org.bson.types.ObjectId;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.Module;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.layers.communication.mind.Mind;
import org.clothocad.core.persistence.Persistor;
import org.mozilla.javascript.NativeArray;
import org.mozilla.javascript.NativeObject;

/**
 *
 * @author spaige
 */
//XXX: actually JavaScriptAPI
public class ScriptAPI {
    ServerSideAPI api;
    Mind mind;
    Persistor persistor;
    
    
    public ScriptAPI(Mind mind, Persistor persistor, String requestId){
        api = new ServerSideAPI(mind, persistor, requestId);
        this.mind = mind;
        this.persistor = persistor;
    }

    public ObjectId create(Map<String, Object> obj) {
        return api.create(new HashMap(obj)); 
    }

    public Map<String, Object> get(Object o) {
        return convertToNative(api.get(o)); 
    }

    public Map<String, Object> get(ObjectId id) {
        return convertToNative(api.get(id)); 
    }

    public List<Map<String, Object>> query(Map<String, Object> spec) {
        return convertToNative(api.query(spec)); 
    }
    
    
    public Object set(Map<String, Object> spec){
        api.set(new HashMap(spec));
        return get(spec.get("id"));
    }
    
    public void destroy(String id){
        api.destroy(id);
    }
    
    public Object run(Object selector, Object arg) throws ScriptException, IllegalAccessException, IllegalArgumentException, InvocationTargetException{
        return run(selector, Arrays.asList(arg));
    }
    
    public Object run(Object selector, List<Object> args) throws ScriptException, IllegalAccessException, IllegalArgumentException, InvocationTargetException{
        Map<String, Object> data = new HashMap<>();
        data.put("id", selector);
        data.put("args", args);
        return api.run(data);
    }
    
    public Object load(Object selector) throws ScriptException {
        boolean isFunction;
        //this, like all the other schema hacks around, needs to be replaced w/ a more sophisticated type-tracking system
        isFunction = api.get(selector).get("schema").toString().endsWith("Function");
        Module module;
        if (isFunction){
            module = persistor.get(Function.class, persistor.resolveSelector(selector.toString(), Function.class, true));
        } else {
            module = persistor.get(Module.class, persistor.resolveSelector(selector.toString(), Module.class, true));
        }
        
        //some kind of circular dependency detection would be good
        return mind.getEngine().eval(module.getCodeToLoad(), module.getLanguage(), this);
    }
    
    
    //XXX: augh, would be best if we had scriptengines that could treat maps as native objects
    //TODO: handle multiple languages
    //XXX: doesn't reach values hidden by non-Map/List fields
    private Map<String,Object> convertToNative(Map<String,Object> obj){
        NativeObject nobj = new NativeObject();
        for (Map.Entry<String, Object> entry : obj.entrySet()) {
            //nobj.defineProperty(entry.getKey(), convertToNative(entry.getValue()), NativeObject.READONLY);
            nobj.put(entry.getKey(), nobj, convertToNative(entry.getValue()));
        }
        
        return nobj;
    }
    
    private Object convertToNative(Object object){
        if (object instanceof Map) return convertToNative((Map) object);
        if (object instanceof List) return convertToNative((List) object);
        return object;
    }
    
    private List convertToNative(List list){
        List convertedObjects = new ArrayList();
        for (Object o : list){
            convertedObjects.add(convertToNative(o));
        }
        
        NativeArray narray = new NativeArray(convertedObjects.toArray());
        return narray;
    }
}

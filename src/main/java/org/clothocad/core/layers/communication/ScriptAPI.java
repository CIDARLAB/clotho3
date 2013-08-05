/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.layers.communication;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.bson.types.ObjectId;
import org.clothocad.core.layers.communication.mind.Mind;
import org.clothocad.core.persistence.Persistor;
import sun.org.mozilla.javascript.Scriptable;
import sun.org.mozilla.javascript.Context;
import sun.org.mozilla.javascript.NativeArray;
import sun.org.mozilla.javascript.NativeObject;

/**
 *
 * @author spaige
 */
public class ScriptAPI {
    ServerSideAPI api;
    
    
    public ScriptAPI(Mind mind, Persistor persistor, String requestId){
        api = new ServerSideAPI(mind, persistor, requestId);
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
    
    public void destroy(ObjectId id){
        api.destroy(id);
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

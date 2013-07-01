/*
Copyright (c) 2010 The Regents of the University of California.
All rights reserved.
Permission is hereby granted, without written agreement and without
license or royalty fees, to use, copy, modify, and distribute this
software and its documentation for any purpose, provided that the above
copyright notice and the following two paragraphs appear in all copies
of this software.

IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.

THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
ENHANCEMENTS, OR MODIFICATIONS.
 */

package org.clothocad.core.persistence;

import java.net.UnknownHostException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import javax.inject.Inject;
import javax.inject.Singleton;
import javax.validation.ConstraintViolation;
import javax.validation.ConstraintViolationException;
import javax.validation.Validation;
import javax.validation.Validator;
import lombok.Delegate;
import lombok.extern.slf4j.Slf4j;
import org.bson.BSONObject;
import org.bson.types.ObjectId;
import org.clothocad.core.aspects.JSONSerializer;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.Sharable;
import org.reflections.Reflections;
import static java.lang.reflect.Modifier.isAbstract;
import java.util.ArrayList;
import java.util.Collection;
import org.clothocad.core.schema.BuiltInSchema;
import org.clothocad.model.Part;

/**
 * @author jcanderson
 * 
 * Manages writing/reading objects to/from the DB, and will eventually also coordinate that with caching
 * Right now is a very simple pass-through class.
 * 
 * queries do not check the cache
 * 
 * fooAsBSON does not check or populate the cache - maintaining a separate BSON cache would be complicating but might help performance-wise
 * 
 * should probably get the cache to cooperate with the object deserializer (or pass the cache in for query methods)
 * possible refactoring: separate deserializer from connection?
 * 
 * having multiple threads have access to the same object sounds like a recipe for disaster -
 * Persistor should hand out copies
 * 
 *   
 */

//TODO: thread safety
//TODO: check out date created/modified/accessed bugs
//TODO: move backend-agnostic logic into persistor
@Singleton
@Slf4j
public class Persistor{
    //TODO: figure out if references to BSON are okay or should be eradicated in favor of Map
    
    private static final int SEARCH_MAX = 5000;
    
    private ClothoConnection connection;
    
    @Delegate
    private JSONSerializer serializer;
    
    private Validator validator = Validation.buildDefaultValidatorFactory().getValidator();
    
    
    @Inject
    public Persistor(final ClothoConnection connection, JSONSerializer serializer){
        this.connection = connection;
        connect();
        this.serializer = serializer;
        
        initializeBuiltInSchemas();
    }
    
    private void validate(ObjBase obj){
        Set<ConstraintViolation<?>> violations = new HashSet<>();
        
        for (ObjBase o : getObjBaseSet(obj)){
            Set<ConstraintViolation<ObjBase>> cvs = validator.validate(o); //XXX: will only validate the current classes
            for (ConstraintViolation violation : cvs){
                violations.add(violation);
            }
        }
        
        if (violations.size() > 0){
            throw new ConstraintViolationException(violations);
        }
    }
    
    
    public <T extends ObjBase> T get(Class<T> type, ObjectId id){
        T obj = connection.get(type, id);
        validate(obj);
        return obj;
    }
    
    //throws ConstraintViolationException, OverwriteConfirmationException
    
    public void save(ObjBase obj) {
        save(obj, false);
    }
    
    public void save(ObjBase obj, boolean overwrite) {
        Set<ObjBase> relevantObjects = getObjBaseSet(obj);
        validate(obj);
        
        if (!overwrite){
            Set<ObjBase> modifiedObjects = new HashSet<>();
            for (ObjBase o : relevantObjects){
                if (modified(o)) modifiedObjects.add(o);
            }
            if (modifiedObjects.size() > 0) throw new OverwriteConfirmationException(modifiedObjects);
        }

        connection.saveAll(relevantObjects);
    }
    
    public ObjectId save(Map<String, Object> data) throws ConstraintViolationException, OverwriteConfirmationException{
        //XXX: convert id field so that mongoDB understands it (hackish)
        if (data.containsKey("id") && !data.containsKey("_id")){
            Object id = data.get("id");
            data.remove("id");
            data.put("_id", id);
        }
        else if (!data.containsKey("_id")){
            Object id = new ObjectId();
            data.put("_id", id);
        }
        
        connection.save(data);
        
        return new ObjectId(data.get("_id").toString());
    }
    
    public void delete(ObjectId id){
        connection.delete(id);
    }
    
    public Map<String, Object> getAsJSON(ObjectId uuid){
        return serializer.toJSON(connection.getAsBSON(uuid).toMap());
    }
    
    
    private Set<ObjBase> getObjBaseSet(ObjBase obj){
        return getObjBaseSet(obj, new HashSet<ObjBase>());
    }
        
    private Set<ObjBase> getObjBaseSet(ObjBase obj, Set<ObjBase> exclude){
        boolean newId = obj.getUUID() == null;
        if(newId) obj.setUUID(ObjectId.get());
        exclude.add(obj);        
        
        //recurse on object's children
        for (ObjBase child: obj.getChildren()) {
            if (!exclude.contains(child) && child != null){
                getObjBaseSet(child, exclude);
            }
        }
        
        return exclude;
    }
    
    public boolean has(ObjectId id){
        return connection.exists(id);
    }
    
    
    private void validateBSON(Map<String, Object> obj) throws ConstraintViolationException {
        //TODO: 
        //get schema set
        //validate for all enforcing schemas
    }
    
    
    private boolean changed(ObjBase obj){
        //TODO
        return true;
    }
    
    private boolean modified(ObjBase obj){
        //TODO
        return false;
    }
    
    public void persistFeature(HashMap<String, Integer> StoreGrams, String feature) {
   /*     System.out.println("Stephanie:  features should be stored in a separate spot than ObjBases.  They aren't UUID-based.  They are are feature-word based");
        try {
            //Convert the StoreGrams to a JSONArray in a JSONObject
            JSONObject bolus = new JSONObject();
            JSONArray data = new JSONArray();
            for(String key : StoreGrams.keySet()) {
                Integer count = StoreGrams.get(key);
                JSONObject item = new JSONObject();
                item.put("command", key);
                item.put("count", count);
                data.put(item);
            }
            bolus.put("data", data);
            bolus.put("feature", feature);

            DBObject bson = ( DBObject ) JSON.parse( bolus.toString() );

            //Query the database for an existing feature entry
            HashMap query = new HashMap();
            query.put("feature", feature);
            BSONObject result = connection.getOneAsBSON(query);
            
            ObjectId oid = null;
            if(result!=null) {
                oid = (ObjectId) result.get("_id");
            } else {
                oid = ObjectId.get();
            }
        
        
            //Install the OID object
            bson.put("_id", oid);
            
            //Save it to the database and return the uuid
            connection.save(new BasicDBObject(bson.toMap()));
            
        } catch (Exception ex) {
            ex.printStackTrace();
        }*/
    }

    public HashMap<String, Integer> loadFeature(String feature) {
      /*  try {
            System.out.println("Stephanie:  feature persistence needs to be separated from ObjBase persistence");

            //Query the database for this feature entry
            HashMap query = new HashMap();
            query.put("feature", feature);
            BSONObject result = connection.getOneAsBSON(query);


            HashMap<String, Integer> out = new HashMap<String, Integer>();
            if(result==null) {
                return out;
            }


            //Transfer the data to a well-typed Map

            //need to re-do the parsing out of the database entry
            String array = result.get("data").toString();
            JSONArray ary = new JSONArray(array);
            
            for(int i=0; i<ary.length(); i++) {
                JSONObject item = ary.getJSONObject(i);
                String command = item.getString("command");
                int count = item.getInt("count");
                out.put(command, count);
            }

            return out;
        } catch(Exception err) {
        }*/
        return null;
    }        

    public Iterable<ObjBase> find(Map<String, Object> query) {
        return find(query, SEARCH_MAX);
    }
    
    public Iterable<ObjBase> find(Map<String, Object> query, int hitmax){
        return connection.get(query);
    }

    public List<Map<String, Object>> findAsBSON(Map<String, Object> spec){
        return findAsBSON(spec, 1000);
    }
    
    public List<Map<String, Object>> findAsBSON(Map<String, Object> spec, int hitmax) {
        List<BSONObject> results = connection.getAsBSON(spec);
        List<Map<String, Object>> out = new ArrayList<>();
        for (BSONObject result : results){
            out.add(serializer.toJSON(result.toMap()));
        }
        return out;
    }
    
    public void connect() {
        try {
            connection.connect();
        } catch (UnknownHostException ex) {
            log.error("Could not connect to database", ex);
        }
        
    }

    public void deleteAll() {
        connection.deleteAll();
        initializeBuiltInSchemas();
   }

    private void initializeBuiltInSchemas() {
        //XXX: just built-in models for now
        Reflections models = new Reflections("org.clothocad.model");

        for (Class<? extends ObjBase> c : models.getSubTypesOf(ObjBase.class)){
            Map<String, Object> query = new HashMap<>();
            query.put("_binaryName", c.getCanonicalName());
            if (!isAbstract(c.getModifiers()) && connection.getAsBSON(query).isEmpty()) {
                save(new BuiltInSchema(c));
            } 
        }
    }

    public <T extends ObjBase> Collection<T> getAll(Class<T> aClass) {
        return connection.getAll(aClass);
    }
}

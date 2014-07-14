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

//TODO: move the translation from clothoisms to mongodb into mongodbconnection

package org.clothocad.core.persistence;

import org.clothocad.core.schema.Converters;
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
import lombok.extern.slf4j.Slf4j;
import org.clothocad.core.datums.ObjBase;
import org.reflections.Reflections;
import java.util.ArrayList;
import java.util.Collection;
import javax.persistence.EntityNotFoundException;
import lombok.Getter;
import org.apache.shiro.SecurityUtils;
import org.apache.shiro.authz.AuthorizationException;
import org.apache.shiro.subject.Subject;
import static org.clothocad.core.ReservedFieldNames.*;
import org.clothocad.core.aspects.Interpreter.GlobalTrie;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.jackson.JSONFilter;
import org.clothocad.core.schema.BuiltInSchema;
import org.clothocad.core.schema.ClothoSchema;
import org.clothocad.core.schema.Converter;
import org.clothocad.core.schema.JavaSchema;
import org.clothocad.core.schema.Schema;
import org.clothocad.core.util.JSON;

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
    public static final int SEARCH_MAX = 5000;
    
    @Getter
    private ClothoConnection connection;
    
    private Validator validator = Validation.buildDefaultValidatorFactory().getValidator();
    
    private Converters converters;
    
    private GlobalTrie globalTrie;
    
    @Inject
    public Persistor(final ClothoConnection connection){
        this(connection, true);
    }
    
    public Persistor(final ClothoConnection connection, boolean initializeBuiltins){
        this.connection = connection;
        connect();
        
       // if (initializeBuiltins) initializeBuiltInSchemas();
        
        converters = new Converters();       
        globalTrie = new GlobalTrie(connection.getCompletionData());
    }
    
    protected void validate(ObjBase obj){
        Set<ConstraintViolation<?>> violations = new HashSet<>();
        
        for (ObjBase o : getObjBaseSet(obj)){
            Set<ConstraintViolation<ObjBase>> cvs = validator.validate(o); //XXX: will only validate the constraints on currently instantiated classes
            for (ConstraintViolation violation : cvs){
                log.info("Constraint violation: {}", violation.getMessage());
                violations.add(violation);
            }
        }
        
        if (violations.size() > 0){
            throw new ConstraintViolationException(violations);
        }
    }
    
    
    public <T extends ObjBase> T get(Class<T> type, ObjectId id) throws EntityNotFoundException {
        try {
            checkPriv(id, "view");
        } catch (AuthorizationException e){
            throw new EntityNotFoundException(id.toString());
        }
        T obj = connection.get(type, id);
        if (obj == null) throw new EntityNotFoundException(id.toString());
        validate(obj);
        return obj;
    }
    
    //throws ConstraintViolationException, OverwriteConfirmationException
    
    public ObjectId save(ObjBase obj) {
        return save(obj, false);
    }
    
    private void checkPriv(ObjectId id, String priviliege) throws AuthorizationException {
        if (id == null) throw new IllegalArgumentException("Null ObjectId");
        Subject currentSubject = SecurityUtils.getSubject();
        if (has(id)) {
            if (!currentSubject.isPermitted("data:"+ priviliege + ":" + id.toString())) {
                log.warn("User {} attempted unauthorized {} on object# {}", currentSubject.getPrincipal(), id);
                throw new AuthorizationException("Not authorized.");
            }
        }
    }
    
    public ObjectId save(ObjBase obj, boolean overwrite) {
        if (obj.getId() != null) checkPriv(obj.getId(), "edit");
        validate(obj);
        
        Set<ObjBase> relevantObjects = getObjBaseSet(obj);

        
        if (!overwrite){
            Set<ObjBase> modifiedObjects = new HashSet<>();
            for (ObjBase o : relevantObjects){
                if (modified(o)) modifiedObjects.add(o);
            }
            if (modifiedObjects.size() > 0) throw new OverwriteConfirmationException(modifiedObjects);
        }

        //recurse in persistor
        connection.saveAll(relevantObjects);
        return obj.getId();
    }
    
    public ObjectId save(Map<String, Object> data) throws ConstraintViolationException, OverwriteConfirmationException{
       if (!SecurityUtils.getSubject().isAuthenticated()){
           throw new AuthorizationException("Anonymous users cannot create or edit objects.");
       }
       if (!data.containsKey(ID)){
            Object id = new ObjectId();
            //Our ObjectId class isn't a BSON datatype, so use string representation
            data.put(ID, id.toString());
       } else {
           checkPriv(new ObjectId(data.get(ID)), "edit");
       }
        connection.save(data);
        
        //Update the GlobalTrie with the new object
        globalTrie.put(data);
        
        return new ObjectId(data.get(ID).toString());
    }
    
    public void delete(ObjectId id) throws AuthorizationException{
        checkPriv(id, "delete");
        connection.delete(id);
    }
    
    public Map<String, Object> getAsJSON(ObjectId uuid){
        return getAsJSON(uuid, null);
    }
    
    public Map<String, Object> getAsJSON(ObjectId id, Set<String> fields){
        try {
            checkPriv(id, "view");
        } catch (AuthorizationException e){
            throw new EntityNotFoundException(id.toString());
        }        
        
        Map<String,Object> result = connection.getAsBSON(id, fields);
        if (result == null) throw new EntityNotFoundException(id.toString());
        return result;
    }
    
    
    //XXX: should return set of possible schemas, not child objects
    // then should not be used as currently is in save
    protected Set<ObjBase> getObjBaseSet(ObjBase obj){
        return getObjBaseSet(obj, new HashSet<ObjBase>());
    }
        
    private Set<ObjBase> getObjBaseSet(ObjBase obj, Set<ObjBase> exclude){
        boolean newId = obj.getId() == null;
        if(newId) obj.setId(new ObjectId());
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
    
    
    public void validateBSON(Map<String, Object> obj) throws ConstraintViolationException {
        //get schema
        if (!obj.containsKey(SCHEMA) || obj.get(SCHEMA) == null){
            throw new IllegalArgumentException("Object does not declare a schema.");
        }
        //validate for all enforcing schemas
        Schema schema = get(Schema.class, new ObjectId(obj.get(SCHEMA)));
        try {
            ObjBase objbase = schema.instantiate(obj);
            validate(objbase);
        } catch (ClassNotFoundException e) {
            throw new IllegalArgumentException("Could not validate schema.", e);
        }
        
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

        query = addSubSchemas(query);
        List<ObjBase> result = connection.get(query, hitmax);
        //TODO: also add converted instances
        
        //filter results for permission
        for (ObjBase obj : result){
            try{
                checkPriv(obj.getId(), "view");
            } catch (AuthorizationException e){
                result.remove(obj);
            }
        }
        return result;
    }

    private List<Map<String,Object>> getConvertedData(Schema originalSchema, Set<String> fields){
        List<Map<String, Object>> results = new ArrayList<>();

        for (Schema schema : converters.getConverterSchemas(originalSchema)){
            Map<String, Object> query = new HashMap<>();
            query.put(SCHEMA, schema.getName());
            List<Map<String,Object>> convertibles = connection.getAsBSON(query);
            Converter converter = converters.getConverter(schema, originalSchema);
            for (Map<String,Object> bson : convertibles){
                Map<String,Object> convertedData = convertAsBSON(bson, schema, converter);
                
                results.add(convertedData);
            }
        }
        
        for (Map<String,Object> result : results){
            filterFields(result, fields);
        }
        
        return results;
    }
    
    public static void filterFields(Map<String,Object> value, Set<String> fields){
        for (String field : value.keySet()){
            if (!fields.contains(field)){
                value.remove(field);
            }
        }
    }
    
    //XXX: mutator - maybe should make copy?
    private Map<String, Object> addSubSchemas(Map<String, Object> query) {
        //XXX: needs to be fixed to handle complex schema queries less clunkily
        if (query.containsKey(SCHEMA)) {
            Object originalSchema = query.get(SCHEMA);
            Map<String, Object> schemaQuery = new HashMap();
            List<String> schemaNames = new ArrayList<>();
            if (originalSchema instanceof Map && ((Map) originalSchema).containsKey("$in")){
                for (Object entry : (List) ((Map) originalSchema).get("$in")){
                    schemaNames.add(entry.toString());
                }
            } else {
                schemaNames.add(originalSchema.toString());
            }
            
            List<String> relatedSchemas = new ArrayList<>();
            
            for (String name : schemaNames) {
                relatedSchemas.addAll(getRelatedSchemas(name));
            }

            if (relatedSchemas.size() > 1) {
                schemaQuery.put("$in", relatedSchemas);
                query.put(SCHEMA, schemaQuery);
            }
        }
        return query;
    }


    //Class set
// an instance's set is any class it has ever been saved as and any classes it was initialized with
//  -> maintaining this set is the connector's responsibility
//if a class T is in an instances set
// all of T's interfaces and T's supertype are in that instance's set
//  -> persistor's responsibility
// (ie, isassignablefrom)

//furthermore, if a converter that accepts T and produces S is available, S is in any set that includes T

//enforcing set is a subset of classes that an instance has been saved as
//forces validation on save IN ADDITION TO the current class the instance is
// -> hm, this could make it impossible to save an instance if the current class and an enforcing class conflict
// -> a good client should have some kind of saveAsNew
//requires write privilege to change (different from instance schema set in that regard)
//defaults to the first class an instance was saved as (?)
    
    public List<Map<String, Object>> findAsBSON(Map<String, Object> spec){
        return findAsBSON(spec, null, 1000);
    }
    
    public List<Map<String, Object>> findAsBSON(Map<String, Object> spec, Set<String> fields, int hitmax) {
        //TODO: limit results
        spec = addSubSchemas(spec);
        List<Map<String,Object>> out = connection.getAsBSON(spec, hitmax, fields);

        
        if (spec.containsKey(SCHEMA) && spec.get(SCHEMA)!= null){
            //copy spec so we don't mutate original search object
            spec = new HashMap(spec);
            Set<String> schemaIds = new HashSet<>();
            Object className = spec.get(SCHEMA).toString();
            //XXX: mongo-specific syntax
            if (className instanceof Map && ((Map) className).containsKey("$in")){
                for (String name : (List<String>) ((Map) className).get("$in")){
                    schemaIds.add(name);
                }
            } else schemaIds.add(className.toString());
            
            for (String id : schemaIds) {
                //try finding a schema by binary name
                Map schemaQuery = new HashMap();
                schemaQuery.put(ID, id);
                Schema originalSchema = connection.getOne(Schema.class, schemaQuery);
                if (originalSchema != null) {
                    List<Map<String, Object>> convertedData = getConvertedData(originalSchema, fields);
                    out.addAll(filterDataByQuery(convertedData, spec));
                }
            }
        }
        out = filterDuplicatesById(out);
        
        //filter results for permission
        for (Map<String,Object> obj : out){
            try{
                checkPriv(new ObjectId(obj.get(ID)), "view");
            } catch (AuthorizationException e){
                out.remove(obj);
            }
        }
        
        return out;
    }
    
    private List<Map<String,Object>> filterDuplicatesById(List<Map<String,Object>> objects){
        List<Map<String,Object>> filteredObjects = new ArrayList<>();
        Set<String> ids = new HashSet<>();
        for (Map<String,Object> object : objects){
            //these objects are coming out of the DB, so we know they have an id
            if (ids.contains(object.get(ID).toString())) continue;
            else {
                ids.add(object.get(ID).toString());
                filteredObjects.add(object);
            }
        }
        return filteredObjects;
    }
   
    private List<Map<String,Object>> filterDataByQuery(List<Map<String,Object>> convertedData,Map<String, Object> spec){
        //TODO
        return convertedData;
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

    protected void initializeBuiltInSchemas() {
        //XXX: just built-in models for now
        Reflections models = new Reflections("org.clothocad");

        for (Class<? extends ObjBase> c : models.getSubTypesOf(ObjBase.class)){
            //XXX: this lets people inject their own implementations of various built-in schemas if they know the name and have db access
            //not sure if this is a problem
            
            if (c.getSuperclass() == ObjBase.class){
                makeBuiltIns(c, null, models);
            }
        }
    }
    
    //XXX: remove this method from 'core' persistor behavior
    private void makeBuiltIns(Class<? extends ObjBase> c, BuiltInSchema superSchema, Reflections ref){
        //JSONfilter only transiently holds data, don't make a schema for it
        if (c == JSONFilter.class) return;
        //XXX: inefficient
        BuiltInSchema builtIn = new BuiltInSchema(c, superSchema);
        if (!has(new ObjectId(c.getCanonicalName()))){
            save(builtIn);           
        }
        for (Class<? extends ObjBase> subClass : ref.getSubTypesOf(c)){
            if (subClass.getSuperclass() == c){
                makeBuiltIns(subClass, builtIn, ref);                
            }
        }
    }
    
    //XXX: doesn't appear to be used anywhere. Remove?
    public <T extends ObjBase> Collection<T> getAll(Class<T> aClass) {
        return connection.getAll(aClass);
    }

    private static final List<Class<? extends Schema>> authoredSchemas = new ArrayList<>();
    static {
        authoredSchemas.add(ClothoSchema.class);
        authoredSchemas.add(JavaSchema.class);
    }
    
    private List<String> getRelatedSchemas(String originalSchemaName) {
        List<String> out = new ArrayList();
        try {
            //get any built-in classes
            //TODO: cache reflections
            Class<?> type = Class.forName(originalSchemaName);
            Reflections r = new Reflections("org.clothocad"); 
            for (Class c : r.getSubTypesOf(type)) {
                out.add(c.getName());
            }
            
            //get any authored schemas
            Schema originalSchema = resolveSchemaFromClassName(originalSchemaName);
            for (Class<? extends Schema> c : authoredSchemas){
                for (Schema schema : this.getAll(c)){
                    if (schema.childOf(originalSchema)){
                        out.add(schema.getBinaryName());
                    }
                }
            }
            

        } catch (ClassNotFoundException | RuntimeException ex) {
        }
        out.add(originalSchemaName);
        return out;
    }

    
    public Map<String,Object> convertAsBSON(Map<String,Object> data, Schema type, Converter converter){
        ObjBase converted = (ObjBase) converter.convert(data, type);
        Map<String,Object> convertedData = JSON.mappify(converted);
        return unifyObjectDescriptions(data, convertedData);
    }
    public static Map<String,Object> unifyObjectDescriptions(Map<String,Object> sourceJSON, Map<String,Object> newJSON) {
        Set<String> exclude = new HashSet();
        exclude.add("schema");
        for (String field : sourceJSON.keySet()){
            if (exclude.contains(field)) continue;
            if (newJSON.containsKey(field)) continue;
            
            newJSON.put(field, sourceJSON.get(field));
        }
        return newJSON;
    }
    
    public Schema resolveSchemaFromClassName(String className){
        Map<String,Object> query = new HashMap();
        query.put("binaryName",className);
        return connection.getOne(Schema.class, query);
    }
    
 
    public List<Map> getCompletions(String word){
        return globalTrie.getCompletions(word);
    }
}

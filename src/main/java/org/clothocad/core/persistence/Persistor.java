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

import com.mongodb.BasicDBObject;
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
import lombok.Delegate;
import lombok.extern.slf4j.Slf4j;
import org.bson.BSONObject;
import org.bson.types.ObjectId;
import org.clothocad.core.aspects.JSONSerializer;
import org.clothocad.core.datums.ObjBase;
import org.reflections.Reflections;
import java.util.ArrayList;
import java.util.Collection;
import javax.persistence.EntityNotFoundException;
import javax.persistence.NonUniqueResultException;
import org.clothocad.core.schema.BuiltInSchema;
import org.clothocad.core.schema.ClothoSchema;
import org.clothocad.core.schema.Converter;
import org.clothocad.core.schema.JavaSchema;
import org.clothocad.core.schema.Schema;
import org.clothocad.model.BasicPartConverter;

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
    
    private Converters converters;
    
    
    @Inject
    public Persistor(final ClothoConnection connection, JSONSerializer serializer){
        this(connection, serializer, true);
    }
    
    public Persistor(final ClothoConnection connection, JSONSerializer serializer, boolean initializeBuiltins){
        this.connection = connection;
        connect();
        this.serializer = serializer;
        
        if (initializeBuiltins) initializeBuiltInSchemas();
        
        
        converters = new Converters();
       
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
    
    
    public <T extends ObjBase> T get(Class<T> type, ObjectId id){
        T obj = connection.get(type, id);
        validate(obj);
        return obj;
    }
    
    //throws ConstraintViolationException, OverwriteConfirmationException
    
    public ObjectId save(ObjBase obj) {
        return save(obj, false);
    }
    
    public ObjectId save(ObjBase obj, boolean overwrite) {
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
        return obj.getUUID();
    }
    
    public ObjectId save(Map<String, Object> data) throws ConstraintViolationException, OverwriteConfirmationException{
        //XXX: convert id field so that mongoDB understands it (hackish)
        if (data.containsKey("id") && !data.containsKey("_id")){
            Object id = data.get("id");
            if (!(id instanceof ObjectId) && ObjectId.isValid(id.toString())){
                id = new ObjectId(id.toString());
            }
            //TODO: figure out what our range of acceptable ids actually is/ what we disable for non-ObjectId ids
            data.remove("id");
            data.put("_id", id);
            
        }
        else if (!data.containsKey("_id")){
            Object id = new ObjectId();
            data.put("_id", id);
        }
        
        //XXX: set up through refs instead
        if (data.containsKey("schema")){
            Object schema = data.get("schema");
            String resolvedSchemaName;
            
            try{
                resolvedSchemaName = get(BuiltInSchema.class,resolveSelector(schema.toString(), true)).getBinaryName();
                if (resolvedSchemaName == null){
                    throw new EntityNotFoundException();
                }
            }
            catch (EntityNotFoundException e){
                try{
                    Map<String,Object> schemaData = getAsJSON(resolveSelector(schema.toString(), false));
                    resolvedSchemaName = schemaData.get("name").toString();
                }
                catch (EntityNotFoundException ex){
                    resolvedSchemaName = schema.toString();
                }
            }
            data.remove("schema");
            data.put("className", resolvedSchemaName);
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
    
    
    //XXX: should return set of possible schemas, not child objects
    // then should not be used as currently is in save
    protected Set<ObjBase> getObjBaseSet(ObjBase obj){
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

        query = modifyQueryForSchemaSearch(query);
        List<ObjBase> result = connection.get(query);
        //TODO: also add converted instances
        return result;
    }

    private List<Map<String,Object>> getConvertedData(Schema originalSchema){
        List<Map<String, Object>> results = new ArrayList<>();

        for (Schema schema : converters.getConverterSchemas(originalSchema)){
            Map<String, Object> query = new HashMap<>();
            query.put("className", schema.getName());
            List<BSONObject> convertibles = connection.getAsBSON(query);
            Converter converter = converters.getConverter(schema, originalSchema);
            for (BSONObject bson : convertibles){
                Map<String,Object> convertedData = convertAsBSON(bson.toMap(), schema, converter);
                
                results.add(convertedData);
            }
        }
        
        return results;
    }
    
    private Map<String, Object> modifyQueryForSchemaSearch(Map<String, Object> query) {
        //if searching for schema
        //XXX: needs to be fixed to handle complex schema queries less clunkily
        if (query.containsKey("className")) {
            Object originalSchema = query.get("className");
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
                query.put("className", schemaQuery);
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
        return findAsBSON(spec, 1000);
    }
    
    public List<Map<String, Object>> findAsBSON(Map<String, Object> spec, int hitmax) {
        spec = modifyQueryForSchemaSearch(spec);
        //TODO: limit results
        List<BSONObject> results = connection.getAsBSON(spec);
        List<Map<String, Object>> out = new ArrayList<>();
        for (BSONObject result : results){
            out.add(serializer.toJSON(result.toMap()));
        }
        
        //also converted stuff
        if (spec.containsKey("className")){
            List<String> schemaNames = new ArrayList<>();
            Object className = spec.get("className");
            if (className instanceof String) schemaNames.add((String) className);
            else if (className instanceof Map && ((Map) className).containsKey("$in")){
                for (String name : (List<String>) ((Map) className).get("$in")){
                    schemaNames.add(name);
                }
            }
            for (String schemaName : schemaNames) {
                //try finding a schema by binary name
                Map schemaQuery = new HashMap();
                schemaQuery.put("binaryName", schemaName);
                Schema originalSchema = connection.getOne(Schema.class, schemaQuery);
                //try finding a schema by name name
                if (originalSchema == null) {
                    originalSchema = get(Schema.class, resolveSelector(schemaName, false));
                }
                List<Map<String, Object>> convertedData = getConvertedData(originalSchema);
                out.addAll(filterDataByQuery(convertedData, spec));
                //TODO: filtering so things are unique
            }
        }
        return filterById(out);
    }
    
    private List<Map<String,Object>> filterById(List<Map<String,Object>> objects){
        List<Map<String,Object>> filteredObjects = new ArrayList<>();
        Set<String> ids = new HashSet<>();
        for (Map<String,Object> object : objects){
            if (ids.contains(object.get("id").toString())) continue;
            else {
                ids.add(object.get("id").toString());
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
    
    private void makeBuiltIns(Class<? extends ObjBase> c, Schema superSchema, Reflections ref){
        Schema builtIn = new BuiltInSchema(c, superSchema);
        try{
            resolveSelector(builtIn.getName(), BuiltInSchema.class, false);
        } catch (EntityNotFoundException e){
            save(builtIn);           
        }
        for (Class<? extends ObjBase> subClass : ref.getSubTypesOf(c)){
            if (subClass.getSuperclass() == c){
                makeBuiltIns(subClass, builtIn, ref);                
            }
        }
    }
    
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
        Map<String,Object> convertedData = serializer.toJSON(converted);
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
        return connection.getOne(Schema.class, new BasicDBObject("binaryName", className));
    }
    
    public ObjectId resolveSelector(String selector, boolean strict) {
        return resolveSelector(selector, null, strict);
    }

    public ObjectId resolveSelector(String selector, Class<? extends ObjBase> type, boolean strict) {
        //uuid?
        ObjectId id;

        try {
            id = new ObjectId(selector);
            return id;
        } catch (IllegalArgumentException e) {
        }

        Map<String, Object> spec = new HashMap<>();
        spec.put("name", selector);
        /* TODO: class & superclass discrimination
         * if (type != null) {
            spec.put("className", type);
        }*/

        //name of something?
        List<Map<String, Object>> results = findAsBSON(spec);

        if (results.isEmpty()) throw new EntityNotFoundException(selector);
        
        id = new ObjectId(results.get(0).get("id").toString());
        if (results.size() == 1 || !strict) {
            return id;
        }

        //complain about ambiguity
        log.warn("Unable to strictly resolve selector {}", selector);
        throw new NonUniqueResultException();
    }
}

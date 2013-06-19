package org.clothocad.core.layers.persistence.mongodb;

import java.lang.reflect.Array;
import java.net.UnknownHostException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.bson.BSONObject;
import org.bson.types.ObjectId;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.util.ClothoDate;
import org.clothocad.core.layers.persistence.ClothoConnection;

import com.github.jmkgreen.morphia.Datastore;
import com.github.jmkgreen.morphia.DatastoreImpl;
import com.github.jmkgreen.morphia.Morphia;
import com.github.jmkgreen.morphia.mapping.DefaultMapper;
import com.github.jmkgreen.morphia.mapping.Mapper;
import com.github.jmkgreen.morphia.mapping.MapperOptions;
import com.google.gson.Gson;
import com.mongodb.BasicDBObject;
import com.mongodb.DB;
import com.mongodb.DBCollection;
import com.mongodb.DBCursor;
import com.mongodb.DBObject;
import com.mongodb.MongoClient;
import com.mongodb.MongoException;
import com.mongodb.WriteResult;
import com.mongodb.util.JSON;
import java.util.Arrays;
import java.util.Date;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.UUID;
import java.util.logging.Level;
import org.bson.BasicBSONObject;
import org.json.JSONException;

import org.json.JSONObject;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 *
 * @author spaige
 */
public class MongoDBConnection
        implements ClothoConnection {
    

    private String host = "localhost";
    private int port = 27017;
    private String dbName = "clotho";
    //the demo will break if you change this without changing the Entity annotation on ObjBase
    private String dataCollName = "data";
    //initialization should be revisited when we integrate parts
    private static Morphia morphia;

    public MongoDBConnection() 
    		throws UnknownHostException {
    	MapperOptions opts = new MapperOptions();
    	opts.objectFactory = new PolymorphicObjectFactory();
    	morphia = new Morphia(new DefaultMapper(opts));
    	
    	this.connect();
    }
    
    private MongoClient connection;
    private DB db;
    private DBCollection data;
    private Datastore dataStore;
    private Mapper mapper;

    @Override
    public void connect() throws UnknownHostException {

        connection = new MongoClient(host, port);
        db = connection.getDB(dbName);
        data = db.getCollection(dataCollName);
        dataStore = new DatastoreImpl(morphia, connection, dbName);
        mapper = dataStore.getMapper();
    }

    @Override
    public boolean isAClothoDatabase() {
        return db != null;
    }

    @Override
    public void disconnect() {
        connection.close();
        connection = null;
        db = null;
        data = null;
        dataStore = null;
    }

    @Override
    public boolean isConnected() {
        return db != null;
    }

//var inst = {"city":"Paris","className":"org.clothocad.model.Institution","country":"United States of America","isDeleted":false,"lastModified":{"$date":"2013-06-08T02:48:10.254Z"},"name":"Test institution","state":"SC"};clotho.say(inst.city);   var obj = clotho.create(inst);
    /**
     * JCA:  This is a hard save to the database with no verification of the data or authorization.  Those checks are done by the calling Aspect.  However, the setting of the UUID is done here.
     * @param json
     * @return H
     */
    public String save(JSONObject json) {
        try {
            System.out.println("JCA implemented, Stephanie change this probably [MongoDBConnection.save] -> "+json);
            
            //Figure out the uuid
            ObjectId id = null;
            if(json.has("id")) {
                id = new ObjectId(json.getString("id"));  //If there is an id field, reuse that UUID
            } else {
                id = ObjectId.get();                      //Otherwise create a new UUID
            }
            
            //Destroy any other UUID references
            json.put("id", "");
            json.put("uuid", "");
            json.put("_id", "");
            json.put("$clotho", "");
            json.remove("id");
            json.remove("uuid");
            json.remove("_id");
            json.remove("$clotho");
            
            //Convert it to a BSON
            String jsonstr = json.toString();
            DBObject bson = ( DBObject ) JSON.parse( jsonstr );

            //Install the OID object
            bson.put("_id", id);
            
            //Save it to the database and return the uuid
            data.save(new BasicDBObject(bson.toMap()));
            return id.toString();
        } catch (Exception ex) {
            ex.printStackTrace();
            return null;
        }
    }
    
    @Override
    //Cascade save
    public void save(ObjBase obj){
        System.out.println("Stephanie:  save in connection should either return false or throw exception when fails.  change required in interface.");
        save(obj, new HashSet<ObjBase>());
    }
    
    //actual meat of cascade save in here
    private boolean save(ObjBase obj, Set<ObjBase> exclude){
    	
        boolean newId = obj.getUUID() == null;
        if(newId) obj.setUUID(ObjectId.get());
        exclude.add(obj);        
        
        //recurse on object's children
        for (ObjBase child: obj.getChildren()) {
            //TODO: needs to detect and skip proxies if lazy loading is used
            if (!exclude.contains(child) && child != null){
                if (!save(child, exclude)){
                    return false;
                }
            }
        }
        
        //write obj to DB
        try {
            singleSave(obj);
        }
        catch (Exception e){
            //if we set an id in anticipation of saving, but the save failed, revert to null (so id is consistent w/ db state)
            if (newId) obj.setUUID(null);
            return false;
        }
        return true;
    }
    
    public void singleSave(ObjBase obj) {
        //needs to update lastUpdated field only if save succeeds
        //needs to check if thing actually needs saving
        //needs to validate object
        obj.setLastModified(new Date());
        if (null != data.findOne(new BasicDBObject("_id", obj.getUUID()))) {
            dataStore.merge(obj);
        } else {
            dataStore.save(obj);
        }

    }

    @Override
    public int save(Collection<ObjBase> objs) {
        int i = 0;
        for (ObjBase o : objs) {
            try{
                save(o);
                i++;
            } catch (Exception e) {
                //logger.error("Error while saving collection", e);
                //aggregate errors and pass back to user - not sure on details yet
            }
        }
        return i;
    }
    

    //non-cascade save
    public void save(BSONObject obj) {
        save(obj.toMap());
    }

    @Override
    public void save(Map obj) {
        data.save(new BasicDBObject(obj));
    }
    
    
    @Override
    //TODO: check for references in database
    //TODO: move to 'deleted' collection instead of hard deleting
    public void delete(ObjBase obj) {
            dataStore.delete(obj);
    }

    @Override
    public int delete(Collection<ObjBase> objs) {
        int i = 0;
        for (ObjBase o : objs) {
            try {
                delete(o);
                i++;
            }
            catch (Exception e) {
                //logger.error("Error while deleting object in collection", e);
                //aggregate errors and pass back to user - not sure on details yet 
            }
        }
        return i;
    }

    @Override
    public Date getTimeModified(ObjBase obj) {
        //TODO: just fetch LastModified field instead of entire object
        ObjBase result = dataStore.get(obj);
        return result.getLastModified();
    }

    @Override
    public <T extends ObjBase> T get(Class<T> type, ObjectId uuid) {
        return dataStore.get(type, uuid);
    }

    @Override
    public BSONObject getAsBSON(ObjectId uuid) {
        return data.findOne(new BasicDBObject("_id", uuid));
    }

    /*Query construction forthcoming
     * @Override
     public ObjBase getOne(BasicDBObject query) {
     Collection<ObjBase> results = get(query);
     for (ObjBase obj : results){
     return obj;
     }
     return null;
     }
     */
    

    public List<ObjBase> get(BSONObject query){
        return get(query.toMap());
    }
    
    @Override
    public List<ObjBase> get(Map query) {
        List<ObjBase> results = new ArrayList<ObjBase>();
        DBCursor cursor = data.find(new BasicDBObject(query));

        try {
            while (cursor.hasNext()) {
                results.add((ObjBase) mapper.fromDBObject(null, cursor.next(), mapper.createEntityCache()));
            }
        } finally {
            cursor.close();
        }
        return results;
    }    

    public List<BSONObject> getAsBSON(BSONObject query) {
        return getAsBSON(query.toMap());
    }

    @Override    
    public List<BSONObject> getAsBSON(Map query) {
        DBCursor cursor = data.find(new BasicDBObject(query));
        List<BSONObject> result = new ArrayList<BSONObject>();
        while (cursor.hasNext()){
            result.add(cursor.next());
        }
        
        return result;
    }
    
    @Override 
    public <T extends ObjBase> List<T> get(Class<T> type, Map query) {
        //TODO:
        return null;
    }

    @Override 
    public <T extends ObjBase> List<T> get(Class<T> type, String name) {
        return dataStore.find(type, "name", name).asList();
    }
    
    @Override    
    public <T extends ObjBase> T getOne(Class<T> type, Map query) {
        DBObject dbResult = data.findOne(new BasicDBObject(query));
        return (T) mapper.fromDBObject(type, dbResult, mapper.createEntityCache());
    }

    @Override
    public BSONObject getOneAsBSON(Map query) {
        return data.findOne(new BasicDBObject(query));
    }
    
    
    @Override
    public void deleteAll() {
        this.data.drop();
    }
    
    public ObjectId[] chain(Map query){
        DBCursor results = data.find(new BasicDBObject(query), new BasicDBObject());
        //TODO: could cause problems with large queries - should chunk w/ skip & limit
        ObjectId[] output = new ObjectId[results.length()];
        for (int i=0;results.hasNext();i++){
            DBObject result = results.next();
            output[i] = (ObjectId) result.get("_id");
        }
        
        return output;
    }
    
    public List<ObjBase> recurseQuery(Map query, String[] path){
        return recurseQuery(query, path, new ArrayList<ObjectId>());
    }
    
    
    //NB: this mutates query
    public List<ObjBase> recurseQuery(Map query, String[] path, List<ObjectId> ids){
        boolean go = true;
        setPath(query, path, ids);
        int origSize;
        while (go){
            origSize = ids.size();
            ids.addAll(Arrays.asList(chain(query)));
            go = ids.size() > origSize;
        }
        
        return get(query);
    }
    
    private static void setPath(Map b, String[] path, Object o){
        if (path.length == 0) return;
        for (int i = 0; i + 1 < path.length; i++){
            b = (Map) b.get(path[i]);
        }
        
        b.put(path[path.length-1], o);
    }

    public void delete(ObjectId id) {
        data.remove(new BasicDBObject("_id", id));           
    }

    public List<ObjBase> get(String name) {
        return get(new BasicDBObject("_id", name).toMap());
    }

    public List<BSONObject> getAsBSON(String name) {
        return getAsBSON(new BasicDBObject("_id", name).toMap());
    }

    public <T extends ObjBase> T getOne(Class<T> type, String name) {
        return getOne(type, new BasicDBObject("_id", name));
    }

    public BSONObject getOneAsBSON(String name) {
        return getOneAsBSON(new BasicDBObject("_id", name));
    }

    @Override
    public <T extends ObjBase> List<T> getAll(Class<T> type) {
        return dataStore.find(type).disableValidation().field(Mapper.CLASS_NAME_FIELDNAME).equal(type.getName()).asList();
        //DBObject dbResult = data.findOne(new BasicDBObject(query.toMap()));
        //return (T) mapper.fromDBObject(type, dbResult, mapper.createEntityCache());
    }

    @Override
    public int saveBSON(Collection<Map> objs) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public <T extends ObjBase> List<BSONObject> getAsBSON(Class<T> type, Map query) {
        query.put(Mapper.CLASS_NAME_FIELDNAME, type.getName());
        return getAsBSON(query);
    }

    @Override
    public <T extends ObjBase> List<BSONObject> getAsBSON(Class<T> type, String name) {
        return getAsBSON(new BasicDBObject("name", name).toMap());
    }

    @Override
    public <T extends ObjBase> BSONObject getOneAsBSON(Class<T> type, Map query) {
        query.put(Mapper.CLASS_NAME_FIELDNAME, type.getName());
        return getOneAsBSON(query);
    }

    @Override
    public <T extends ObjBase> BSONObject getOneAsBSON(Class<T> type, String name) {
        return getOneAsBSON(new BasicDBObject("name", name));
    }

}

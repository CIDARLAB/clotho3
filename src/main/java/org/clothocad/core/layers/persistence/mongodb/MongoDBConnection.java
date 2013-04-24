package org.clothocad.core.layers.persistence.mongodb;

import java.lang.reflect.Array;
import java.net.UnknownHostException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.bson.BSONObject;
import org.bson.types.ObjectId;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.layers.persistence.ClothoConnection;

import com.github.jmkgreen.morphia.Datastore;
import com.github.jmkgreen.morphia.DatastoreImpl;
import com.github.jmkgreen.morphia.Morphia;
import com.github.jmkgreen.morphia.mapping.DefaultMapper;
import com.github.jmkgreen.morphia.mapping.Mapper;
import com.github.jmkgreen.morphia.mapping.MapperOptions;
import com.mongodb.BasicDBObject;
import com.mongodb.DB;
import com.mongodb.DBCollection;
import com.mongodb.DBCursor;
import com.mongodb.DBObject;
import com.mongodb.MongoClient;
import com.mongodb.MongoException;
import java.util.Arrays;
import java.util.Date;
import java.util.HashSet;
import java.util.Set;

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

    {
        {
            MapperOptions opts = new MapperOptions();
            opts.objectFactory = new PolymorphicObjectFactory();
            morphia = new Morphia(new DefaultMapper(opts));
        }
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

    
    
    
    /*
    private boolean save(Set<ObjBase> exclude) {
        boolean newId = this.UUID == null;
        if(newId) this.UUID = ObjectId.get();
        exclude.add(this);        
        List<ObjBase> children = getChildren();
        
        for (ObjBase child: children) {
            //TODO: needs to detect and skip proxies if lazy loading is used
            if (!exclude.contains(child) && child != null){
                if (!child.save(exclude)){
                    return false;
                }
            }
        }
        
        //save self
        
        if  (!Persistor.get().persist(this)) {
            //if we set an id in anticipation of saving, but the save failed, revert to null (so id is consistent w/ db state)
            if (newId) this.UUID = null;
            return false;
        }
        return true;
    }     */
    
    
    
    @Override
    //Cascade save
    public void save(ObjBase obj){
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
                //logger.error();
                //aggregate errors and pass back to user - not sure on details yet
            }
        }
        return i;
    }
    
    @Override
    //non-cascade save
    public void save(BSONObject obj) {
        data.save(new BasicDBObject(obj.toMap()));
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
                //logger.error();
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
    public <T> T get(Class<T> type, ObjectId uuid) {
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
    
    @Override
    public List<ObjBase> get(BSONObject query) {
        List<ObjBase> results = new ArrayList<ObjBase>();
        DBCursor cursor = data.find(new BasicDBObject(query.toMap()));

        try {
            while (cursor.hasNext()) {
                results.add((ObjBase) mapper.fromDBObject(null, cursor.next(), mapper.createEntityCache()));
            }
        } finally {
            cursor.close();
        }
        return results;
    }

    
    @Override    
    public Collection<BSONObject> getAsBSON(BSONObject query) {
        DBCursor cursor = data.find(new BasicDBObject(query.toMap()));
        Collection<BSONObject> result = new ArrayList<BSONObject>();
        while (cursor.hasNext()){
            result.add(cursor.next());
        }
        
        return result;
    }

    @Override    
    public <T> T getOne(Class<T> type, BSONObject query) {
        DBObject dbResult = data.findOne(new BasicDBObject(query.toMap()));
        return (T) mapper.fromDBObject(type, dbResult, mapper.createEntityCache());
    }

    @Override
    public BSONObject getOneAsBSON(BSONObject query) {
        return data.findOne(new BasicDBObject(query.toMap()));
    }
    
    
    @Override
    public void deleteAll() {
        this.data.drop();
    }
    
    public ObjectId[] chain(BSONObject query){
        DBCursor results = data.find(new BasicDBObject(query.toMap()), new BasicDBObject());
        //TODO: could cause problems with large queries - should chunk w/ skip & limit
        ObjectId[] output = new ObjectId[results.length()];
        for (int i=0;results.hasNext();i++){
            DBObject result = results.next();
            output[i] = (ObjectId) result.get("_id");
        }
        
        return output;
    }
    
    public List<ObjBase> recurseQuery(BSONObject query, String[] path){
        return recurseQuery(query, path, new ArrayList<ObjectId>());
    }
    
    
    //NB: this mutates query
    public List<ObjBase> recurseQuery(BSONObject query, String[] path, List<ObjectId> ids){
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
    
    private static void setPath(BSONObject b, String[] path, Object o){
        if (path.length == 0) return;
        for (int i = 0; i + 1 < path.length; i++){
            b = (BSONObject) b.get(path[i]);
        }
        
        b.put(path[path.length-1], o);
    }

}

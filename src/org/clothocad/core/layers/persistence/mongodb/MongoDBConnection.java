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
import com.mongodb.BasicDBObject;
import com.mongodb.DB;
import com.mongodb.DBCollection;
import com.mongodb.DBCursor;
import com.mongodb.DBObject;
import com.mongodb.MongoClient;
import com.mongodb.MongoException;

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
    {{
        MapperOptions opts = new MapperOptions();
        opts.objectFactory = new PolymorphicObjectFactory();
        morphia = new Morphia(new DefaultMapper(opts));
    }}
    
    private MongoClient connection;
    private DB db;
    private DBCollection data;
    private Datastore dataStore;
    private Mapper mapper;
    
    @Override
    public boolean connect() {
        try {
            connection = new MongoClient(host, port);
            db = connection.getDB(dbName);
            data = db.getCollection(dataCollName);
            dataStore = new DatastoreImpl(morphia, connection, dbName);
            mapper = dataStore.getMapper();
        }
        catch (UnknownHostException e){
            return false;
        }
        catch (MongoException e){
            return false;
        }
        return true;
        
    }

    public void clearDB() {
    	
    }
    
    @Override
    public boolean isAClothoDatabase() {
        return db != null;
    }

    @Override
    public boolean disconnect() {
        connection.close();
        connection = null;
        db = null;
        data = null;
        dataStore = null;
        return true;
    }

    @Override
    public boolean isConnected() {
        return db != null;
    }

    @Override
    public boolean save(ObjBase obj) {
        //needs to update lastUpdated field only if save succeeds
        //needs to check if thing actually needs saving
        //needs to validate object
        obj.setLastModified(new ClothoDate());
        if (null != data.findOne(new BasicDBObject("_id", obj.getUUID()))){
            dataStore.merge(obj);
        } else {
            dataStore.save(obj);
        }
            
        return true;
    }

    @Override
    public int save(Collection<ObjBase> objs) {
        int i = 0;
        boolean saved;
        for (ObjBase o : objs){
            saved = save(o);
            if (saved) i ++;
        }
        return i;
    }
    
    public boolean save(BSONObject obj){
        data.save(new BasicDBObject(obj.toMap()));
        return true;
    }
    

    @Override
    //TODO: check for references in database
    //TODO: move to 'deleted' collection instead of just deleting
    public boolean delete(ObjBase obj) {
        try{
            dataStore.delete(obj);
        }
        catch(MongoException e){
            return false;
        }
        return true;
    }

    @Override
    public int delete(Collection<ObjBase> objs) {
        int i = 0;
        boolean deleted;
        for (ObjBase o : objs){
            deleted = delete(o);
            if (deleted) i ++;
        }
        return i;
    }

    @Override
    public ClothoDate getTimeModified(ObjBase obj) {
        //TODO: just fetch LastModified field instead of entire object
        ObjBase result = dataStore.get(obj);
        return result.getLastModified();
    }

    @Override
    public <T> T get(Class<T> type, ObjectId uuid) {
        return dataStore.get(type, uuid);
    }
    
    @Override
    public BSONObject getAsBSON(ObjectId uuid){
        return data.findOne(new BasicDBObject("_id",uuid));
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
             
        try{
            while(cursor.hasNext()){
                results.add((ObjBase) mapper.fromDBObject(null, cursor.next(), mapper.createEntityCache()));
            }
        }
        finally {
            cursor.close();
        }
        return results;
    }
    
    public BSONObject getAsBSON(BSONObject query) {
        DBCursor cursor  = data.find(new BasicDBObject(query.toMap()));
        return new BasicDBObject("results", cursor.toArray());        
    }

    @Override
    public ClothoQuery createQuery(Class type) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public <T> T getOne(Class<T> type, BSONObject query) {
        DBObject dbResult = data.findOne(new BasicDBObject(query.toMap()));
        return (T) mapper.fromDBObject(type, dbResult, mapper.createEntityCache());
    }
    
    
    public ObjBase[] getTableAsArray(ObjBase obj) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public <T> T[] getTableAsArray(Class<T> type) {
		BasicDBObject query = new BasicDBObject("className", type.getCanonicalName());		
		List<ObjBase> lst = this.get(query);

		T[] array = (T[])Array.newInstance(type, lst.size());
		return lst.toArray(array);
    }

	@Override
	public boolean clear() {
		this.data.drop();
		return true;
	}
}

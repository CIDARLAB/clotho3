/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jongo;

import org.clothocad.core.persistence.jackson.JSONFilter;
import static com.fasterxml.jackson.annotation.JsonInclude.Include.NON_NULL;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.core.type.TypeReference;
import static com.fasterxml.jackson.databind.DeserializationFeature.FAIL_ON_UNKNOWN_PROPERTIES;
import com.fasterxml.jackson.databind.ObjectMapper;
import static com.fasterxml.jackson.databind.SerializationFeature.FAIL_ON_EMPTY_BEANS;
import com.mongodb.BasicDBObject;
import com.mongodb.DB;
import com.mongodb.DBCollection;
import com.mongodb.DBObject;
import com.mongodb.MongoClient;
import static com.mongodb.MongoException.DuplicateKey;
import java.net.UnknownHostException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import javax.inject.Inject;
import javax.inject.Named;
import lombok.extern.slf4j.Slf4j;
import org.apache.shiro.authc.SimpleAccount;
import org.apache.shiro.crypto.hash.SimpleHash;
import org.apache.shiro.util.ByteSource;
import static org.clothocad.core.ReservedFieldNames.*;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.ClothoConnection;
import org.clothocad.core.security.CredentialStore;
import org.jongo.ResultHandler;
import org.python.google.common.collect.Lists;

/**
 *
 * @author spaige
 */
@Slf4j
public class JongoConnection implements ClothoConnection, CredentialStore {

    protected RefJongo jongo;
    protected DB db;
    protected MongoClient client;
    protected ObjectMapper mapper;
    //TODO: split into separate collections per top-level schema
    protected RefMongoCollection data;
    protected DBCollection cred;
    protected DBCollection rawDataCollection;
    private static final TypeReference<Map<String, Object>> STRINGMAP = new TypeReference<Map<String, Object>>() {
    };

    @Inject
    public JongoConnection(@Named("dbport") int port, @Named("dbhost") String host, @Named("dbname") String dbname) throws UnknownHostException {
        db = new MongoClient(host, port).getDB(dbname);
        rawDataCollection = db.getCollection("data");
        cred = db.getCollection("cred");

    }

    @Override
    public void connect() throws UnknownHostException {
        //TODO: cover reconnect case?

        //Mimic Jongo customization         
        mapper = new ObjectMapper();
        mapper.disable(FAIL_ON_UNKNOWN_PROPERTIES);
        mapper.setSerializationInclusion(NON_NULL);
        mapper.disable(FAIL_ON_EMPTY_BEANS);
        //mapper.enableDefaultTyping(ObjectMapper.DefaultTyping.OBJECT_AND_NON_CONCRETE);
        //jongo mimicking over
        
        /**
         * redundant with ObjBase annotations
         *
         * mapper.disable(AUTO_DETECT_GETTERS);
         * mapper.disable(AUTO_DETECT_IS_GETTERS);
         **/

        jongo = new RefJongo(db, new ClothoMapper());
        data = jongo.getCollection("data");
    }

    //Do we really need this?
    @Override
    public boolean isAClothoDatabase() {
        //TODO
        return db != null;
    }

    @Override
    public void disconnect() {
        client.close();
    }

    @Override
    public boolean isConnected() {
        return db != null;
    }

    @Override
    public void save(ObjBase obj) {
        try{
            data.save(obj);
        }
        catch (DuplicateKey e){
            data.update("{_id:#}", obj.getId().toString()).with(obj);
        }
    }

    @Override
    public void save(Map obj) {
        obj = mongifyIdField(obj);
        rawDataCollection.save(new BasicDBObject(obj));
    }

    @Override
    public void saveAll(Iterable<ObjBase> objs) {
        for (ObjBase o : objs) {
            save(o);
        }
    }

    @Override
    public int saveBSON(Collection<Map> objs) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void delete(ObjBase obj) {
        //TODO: check for references in database
        //TODO: move to 'deleted' collection instead of hard deleting
        delete(obj.getId());
    }

    @Override
    public void delete(ObjectId id) {
        data.remove("{_id:#}", id.toString());
    }

    @Override
    public int delete(Collection<ObjBase> objs) {
        int i = 0;
        for (ObjBase o : objs) {
            try {
                delete(o);
                i++;
            } catch (Exception e) {
                log.error("Error while deleting object in collection", e);
                //aggregate errors and pass back to user - not sure on details yet 
            }
        }
        return i;
    }

    @Override
    public Date getTimeModified(ObjBase obj) {
        //TODO: just fetch LastModified field instead of /entire object
        // use data.projection('{lastModified:1}').as(...?
        ObjBase result = data.resolvingFindOne("{_id:#}", obj.getId().toString()).as(ObjBase.class);
        return result.getLastModified();
    }

    @Override
    public <T extends ObjBase> T get(Class<T> type, ObjectId uuid) {
        return data.resolvingFindOne("{_id:#}", uuid.toString()).as(type);
    }

    @Override
    public Map<String, Object> getAsBSON(ObjectId uuid) {
        return data.findOne("{_id:#}", uuid.toString()).map(DemongifyHandler.get());
    }

    @Override
    public List<ObjBase> get(Map query) {
        return Lists.newArrayList(data.resolvingFind(serialize(mongifyIdField(query))).as(ObjBase.class));
    }

    @Override
    public List<ObjBase> get(String name) {
        return Lists.newArrayList(data.resolvingFind("{name:#}", name).as(ObjBase.class));
    }

    @Override
    public <T extends ObjBase> List<T> get(Class<T> type, Map query) {
        return Lists.newArrayList(data.resolvingFind(serialize(mongifyIdField(query))).as(type));
    }

    @Override
    public <T extends ObjBase> List<T> get(Class<T> type, String name) {
        return Lists.newArrayList(data.resolvingFind("{name:#}", name).as(type));

    }

    @Override
    public List<Map<String, Object>> getAsBSON(Map query) {
        return Lists.newArrayList(data.find(serialize(mongifyIdField(query))).map(DemongifyHandler.get()));
    }

    @Override
    public List<Map<String, Object>> getAsBSON(String name) {
        return getAsBSON(new BasicDBObject("name", name).toMap());
    }

    @Override
    public <T extends ObjBase> List<Map<String, Object>> getAsBSON(Class<T> type, Map query) {
        throw new UnsupportedOperationException("Not supported yet."); 
    }

    @Override
    public <T extends ObjBase> List<Map<String, Object>> getAsBSON(Class<T> type, String name) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public <T extends ObjBase> T getOne(Class<T> type, Map query) {
        return data.findOne(serialize(mongifyIdField(query))).as(type);
    }

    @Override
    public <T extends ObjBase> T getOne(Class<T> type, String name) {
        return data.findOne("{name:#}", name).as(type);
    }

    @Override
    public Map<String, Object> getOneAsBSON(Map query) {
        return data.findOne(serialize(mongifyIdField(query))).map(DemongifyHandler.get());
    }

    @Override
    public Map<String, Object> getOneAsBSON(String name) {
        return getOneAsBSON(new BasicDBObject("name", name).toMap());
    }

    @Override
    public <T extends ObjBase> Map<String, Object> getOneAsBSON(Class<T> type, Map query) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public <T extends ObjBase> Map<String, Object> getOneAsBSON(Class<T> type, String name) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public <T extends ObjBase> List<T> getAll(Class<T> type) {
        return Lists.newArrayList(data.resolvingFind("{schema:#}", type.getCanonicalName()).as(type));
    }

    @Override
    public void deleteAll() {
        rawDataCollection.drop();
    }

    @Override
    public boolean exists(ObjectId id) {
        return rawDataCollection.find(new BasicDBObject("_id", id.toString())).hasNext();
    }

    @Override
    public SimpleAccount getAccount(String username) {
        DBObject accountData = cred.findOne(new BasicDBObject("_id", username));
        if (accountData == null) {
            return null;
        }

        SimpleAccount account = new SimpleAccount(username, accountData.get("hash"), ByteSource.Util.bytes(accountData.get("salt")), "clotho");
        return account;
    }

    @Override
    public void saveAccount(String username, SimpleHash hashedPw, ByteSource salt) {
        //create account needs to fail if username exists
        DBObject account = new BasicDBObject("_id", username);
        account.put("hash", hashedPw.getBytes());
        account.put("salt", salt.getBytes());

        cred.save(account);
    }

    @Override
    public void deleteAllCredentials() {
        cred.drop();
    }

    private List<Map<String, Object>> mappify(Iterable<JSONFilter> objs) {
        List<Map<String, Object>> out = new ArrayList<>();
        for (JSONFilter obj : objs) {
            out.add(mappify(obj));
        }
        return out;
    }

    private Map<String, Object> mappify(JSONFilter obj) {

        return mapper.convertValue(obj, STRINGMAP);
    }

    private String serialize(Object obj) {
        try {
            return mapper.writeValueAsString(obj);
        } catch (JsonProcessingException ex) {
            throw new RuntimeException(ex);
        }
    }

    //renames "id" key to "_id", and replaces ObjectId with String
    protected static Map<String,Object> mongifyIdField(Map<String, Object> obj) {
        if (obj.containsKey(ID)){
            Map<String,Object> copy = new HashMap<>();
            copy.putAll(obj);
            obj = copy;
            Object id = obj.get(ID);
            if (id instanceof ObjectId){
                id = ((ObjectId) id).toString();
            }
            obj.remove(ID);
            obj.put("_id", id);
        }
        return obj;
    }
    
    //renames "_id" to "id" 
    //doesn't recurse on sub-objects, since no embedded object will have an objectid
    //mutates instead of returning new object
    protected static Map<String,Object> demongifyIdField(Map<String,Object> obj){
        if (obj.containsKey("_id")){
            Object id = obj.get("_id");
            obj.remove("_id");
            obj.put(ID, id);
        }
        return obj;
    }
    
    protected static class DemongifyHandler implements ResultHandler<Map<String,Object>>{

        @Override
        public Map<String,Object> map(DBObject result) {
            BasicDBObject basicResult = new BasicDBObject();
            basicResult.putAll(result);
            Map<String,Object> resultMap = basicResult.toMap();
            return demongifyIdField(resultMap);
        }
        
        private static DemongifyHandler instance;
        
        public static DemongifyHandler get(){
            if (instance == null){
                instance = new DemongifyHandler();
            }
            return instance;
        }

    } 
}

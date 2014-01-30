/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jongo;

import com.mongodb.BasicDBObject;
import com.mongodb.DB;
import com.mongodb.DBCollection;
import com.mongodb.DBCursor;
import com.mongodb.DBObject;
import com.mongodb.MongoClient;
import java.net.UnknownHostException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.List;
import java.util.Map;
import javax.inject.Inject;
import javax.inject.Named;
import lombok.extern.slf4j.Slf4j;
import org.apache.shiro.authc.SimpleAccount;
import org.apache.shiro.crypto.hash.SimpleHash;
import org.apache.shiro.util.ByteSource;
import org.bson.BSONObject;
import org.bson.types.ObjectId;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.persistence.ClothoConnection;
import org.clothocad.core.security.CredentialStore;
import org.jongo.Jongo;
import org.jongo.MongoCollection;
import org.python.google.common.collect.Lists;

/**
 *
 * @author spaige
 */

@Slf4j
public class JongoConnection implements ClothoConnection, CredentialStore {

    protected Jongo jongo;
    protected DB db;
    protected MongoClient client;
    //TODO: split into separate collections per top-level schema
    protected MongoCollection data;
    protected DBCollection cred;
    protected DBCollection rawDataCollection;
    /*
     * Jackson config:
     * serialize all fields (except transient) s'thing like:
     * @JsonAutoDetect(fieldVisibility=Visibility.ANY, getterVisibility=Visibility.NONE, isGetterVisibility=Visibility.NONE, setterVisibility=Visibility.NONE)
     */
    
    @Inject
    public JongoConnection(@Named("dbport") int port, @Named("dbhost") String host, @Named("dbname") String dbname) throws UnknownHostException{
        db = new MongoClient(host, port).getDB(dbname);

    }
    
    @Override
    public void connect() throws UnknownHostException {
    //TODO: cover reconnect case?
        jongo = new Jongo(db);    
    }

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
        data.save(obj);
    }

    @Override
    public void save(Map obj) {
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
        data.remove(id);
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
                log.error("Error while deleting object in collection", e);
                //aggregate errors and pass back to user - not sure on details yet 
            }
        }
        return i;    }

    @Override
    public Date getTimeModified(ObjBase obj) {
        //TODO: just fetch LastModified field instead of /entire object
        // use data.projection('{lastModified:1}').as(...?
        ObjBase result = data.findOne(obj.getId()).as(ObjBase.class);
        return result.getLastModified();
    }

    @Override
    public <T extends ObjBase> T get(Class<T> type, ObjectId uuid) {
        return data.findOne(uuid).as(type);
    }

    @Override
    public BSONObject getAsBSON(ObjectId uuid) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public List<ObjBase> get(Map query) {
        return Lists.newArrayList(data.find(query.toString()).as(ObjBase.class));
    }

    @Override
    public List<ObjBase> get(String name) {
        return Lists.newArrayList(data.find("{name:#}", name).as(ObjBase.class));
    }

    @Override
    public <T extends ObjBase> List<T> get(Class<T> type, Map query) {
        return Lists.newArrayList(data.find(query.toString()).as(type));
    }

    @Override
    public <T extends ObjBase> List<T> get(Class<T> type, String name) {
        return Lists.newArrayList(data.find("{name:#}").as(type));
  
    }

    @Override
    public List<BSONObject> getAsBSON(Map query) {
        DBCursor cursor = rawDataCollection.find(new BasicDBObject(query));
        List<BSONObject> result = new ArrayList<>();
        while (cursor.hasNext()){
            result.add(cursor.next());
        }
        
        return result;    
    }

    @Override
    public List<BSONObject> getAsBSON(String name) {
        return getAsBSON(new BasicDBObject("name", name).toMap());
    }

    @Override
    public <T extends ObjBase> List<BSONObject> getAsBSON(Class<T> type, Map query) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public <T extends ObjBase> List<BSONObject> getAsBSON(Class<T> type, String name) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public <T extends ObjBase> T getOne(Class<T> type, Map query) {
        return data.findOne(query.toString()).as(type);
    }

    @Override
    public <T extends ObjBase> T getOne(Class<T> type, String name) {
        return data.findOne("{name:#}", name).as(type);
    }

    @Override
    public BSONObject getOneAsBSON(Map query) {
        return rawDataCollection.findOne(new BasicDBObject(query));
    }

    @Override
    public BSONObject getOneAsBSON(String name) {
        return rawDataCollection.findOne(new BasicDBObject("name", name));
    }

    @Override
    public <T extends ObjBase> BSONObject getOneAsBSON(Class<T> type, Map query) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public <T extends ObjBase> BSONObject getOneAsBSON(Class<T> type, String name) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public <T extends ObjBase> List<T> getAll(Class<T> type) {
        return Lists.newArrayList(data.find("{schema:#}", type.toString()).as(type));
    }

    @Override
    public void deleteAll() {
        rawDataCollection.drop();
    }

    @Override
    public boolean exists(ObjectId id) {
        return null != rawDataCollection.find(new BasicDBObject("_id", id));
    }

    @Override
    public SimpleAccount getAccount(String username) {
        DBObject accountData = cred.findOne(new BasicDBObject("_id", username));
        if (accountData == null) return null;
        
        SimpleAccount account = new SimpleAccount(username, accountData.get("hash"), ByteSource.Util.bytes(accountData.get("salt")), "clotho");
        return account; 
    }

    @Override
    public void saveAccount(String username, SimpleHash hashedPw, ByteSource salt) {
        //create account needs to fail if username exists
        DBObject account = new BasicDBObject("_id", username);
        account.put("hash", hashedPw.getBytes());
        account.put("salt", salt.getBytes());
        
        cred.save(account);    }

    @Override
    public void deleteAllCredentials() {
        cred.drop();
    }
    
}

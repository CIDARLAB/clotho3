/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.communication;

import com.fasterxml.jackson.core.JsonParseException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import static org.clothocad.core.ReservedFieldNames.*;
import org.clothocad.core.datums.ObjectId;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.Ignore;

/**
 *
 * @author spaige
 */
public class ServerAPITest extends AbstractServerAPITest {


    public ServerAPITest() {
        super();
    }

    @Test
    public void get() {
        Map<String, Object> result = api.get(ids.get(2));
        new ObjectId(((List)result.get("composition")).get(0).toString());

    }

    @Test
    public void getNonExistent() {
        assertEquals(null, api.get(new ObjectId()));
    }

    @Ignore("authz not ready")
    public void getPrivate() {
    }

    @Test
    public void set() {
        //TODO!
    }

    @Test
    public void setNonExistent() {
        //XXX: this is actually bad data
        Map<String, Object> newPart = new HashMap<>();
        newPart.put("name", "Test Part 5");
        newPart.put("description", "previously unsaved test part");
        newPart.put("sequence", "ATCGATCG");
        newPart.put("format", "FreeForm");

        
        
        ObjectId id = api.set(new HashMap(newPart));
        
        Map<String, Object> createdPart = api.get(id);
        
        for (String key : newPart.keySet()){
            assertEquals(newPart.get(key), createdPart.get(key));
        }
    }

    @Test
    public void setInvalid() {
    }

    @Ignore("authz not ready")
    public void setPrivate() {
    }

    @Test
    public void create() {
        ObjectId id = new ObjectId();

        //XXX: this is actually bad data
        Map<String, Object> newPart = new HashMap<>();
        newPart.put("name", "Test Part 4");
        newPart.put("description", "previously unsaved test part");
        newPart.put("sequence", "CCCCCCCCCCCCCCCCCCCCCCCC");
        newPart.put("format", "FreeForm");
        newPart.put(ID, id.toString());

        ObjectId createdId = api.create(newPart);

        assertEquals(id, createdId);
        assertNotNull(api.get(id));
    }

    @Test
    public void createExisting() {
        TestConnection connection = new TestConnection("test");
        mind.setConnection(connection);

        ObjectId id = new ObjectId();

        Map<String, Object> obj = new HashMap();
        obj.put(ID, id);

        api.create(obj);
        api.create(obj);
        
        assertEquals(ServerSideAPI.Severity.FAILURE, ((Map) connection.messages.get(1).getData()).get("class"));

        obj = new HashMap();
        obj.put(ID, id);
        
        api.create(obj);
        
        assertEquals(ServerSideAPI.Severity.FAILURE, ((Map) connection.messages.get(2).getData()).get("class"));
    }

    public void createWithBadId(){
        
    }
    
    @Ignore("authz not ready")
    public void createWithoutPrivs() {
    }
    
    

    @Test
    public void destroy() {
    }

    @Ignore("authz not ready")
    public void destroyWithoutPrivs() {
    }

    @Test
    public void query() throws JsonParseException {
        //TODO: implement schema set
        TestConnection connection = new TestConnection("test");
        mind.setConnection(connection);

        //filter out unseen results
        Map<String, Object> query = new HashMap<>();
        query.put("schema", "org.clothocad.model.Part");
        
        List<Map<String, Object>> results = api.query(query);
        Set<String> names = new HashSet();

        for (Map<String, Object> result : results) {
            names.add(result.get("name").toString());
        }

        assertTrue(names.contains("Test Part 1"));
        assertTrue(names.contains("Test Part 2"));
        assertTrue(names.contains("Test Part 3"));
        
        //assert that sequence is either in the database or embedded
    }

    @Test
    public void run() {
    }

    @Test
    public void runNonExistent() {
    }

    @Test
    public void runWrongArguments() {
    }

    @Test
    public void runExecutionError() {
    }

    @Ignore("authz not ready")
    public void runWithoutPrivs() {
    }
}

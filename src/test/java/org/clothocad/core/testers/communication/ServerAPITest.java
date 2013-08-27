/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.testers.communication;

import com.fasterxml.jackson.core.JsonParseException;
import com.google.inject.Guice;
import com.google.inject.Injector;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.bson.types.ObjectId;
import org.clothocad.core.layers.communication.ServerSideAPI;
import org.clothocad.core.layers.communication.mind.Mind;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.persistence.mongodb.MongoDBModule;
import org.clothocad.core.testers.ClothoTestModule;
import org.clothocad.core.utils.TestUtils;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.Ignore;

/**
 *
 * @author spaige
 */
public class ServerAPITest {

    private static ServerSideAPI api;
    private static ServerSideAPI unprivilegedUser;
    private static Persistor persistor;
    private static List<ObjectId> ids;
    private static Mind mind;

    public ServerAPITest() {
    }

    @BeforeClass
    public static void setUpClass() {
        Injector injector = Guice.createInjector(new ClothoTestModule(), new MongoDBModule());
        persistor = injector.getInstance(Persistor.class);
        mind = new Mind();
        api = new ServerSideAPI(mind, persistor, null);
        persistor.connect();
        mind.setConnection(new TestConnection("test"));
    }

    @AfterClass
    public static void tearDownClass() {
    }

    @Before
    public void setUp() {
        persistor.deleteAll();
        ids = TestUtils.setupTestData(persistor);
    }

    @After
    public void tearDown() {
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
        newPart.put("id", id.toString());

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
        obj.put("id", id);

        api.create(obj);
        api.create(obj);
        
        assertEquals(ServerSideAPI.Severity.FAILURE, ((Map) connection.messages.get(1).data).get("class"));

        obj = new HashMap();
        obj.put("_id", id);
        
        api.create(obj);
        
        assertEquals(ServerSideAPI.Severity.FAILURE, ((Map) connection.messages.get(2).data).get("class"));
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
        //TODO: switch back to Part, implement schema set
        TestConnection connection = new TestConnection("test");
        mind.setConnection(connection);

        //filter out unseen results
        Map<String, Object> query = new HashMap<>();
        query.put("schema", "Part");
        
        List<Map<String, Object>> results = api.query(query);
        assertEquals(4, results.size());
        //assertEquals(3, results.size());
        Set<String> names = new HashSet();

        for (Map<String, Object> result : results) {
            names.add(result.get("name").toString());
        }

        assertTrue(names.contains("Test Part 1"));
        assertTrue(names.contains("Test Part 2"));
        assertTrue(names.contains("Test Part 3"));
        assertTrue(names.contains("B0015"));
        
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
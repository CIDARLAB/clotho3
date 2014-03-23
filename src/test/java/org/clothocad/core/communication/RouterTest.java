/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.communication;

import com.google.inject.Guice;
import com.google.inject.Injector;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.persistence.EntityNotFoundException;
import org.apache.shiro.SecurityUtils;
import org.apache.shiro.mgt.SecurityManager;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.persistence.jongo.JongoModule;
import org.clothocad.core.security.ClothoRealm;
import org.clothocad.core.testers.ClothoTestModule;
import org.clothocad.core.util.JSON;
import org.clothocad.core.util.TestUtils;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author spaige
 */
public class RouterTest {

    private static Router router;
    private static Injector injector;
    private static List<ObjectId> ids;

    public RouterTest() {
    }

    @BeforeClass
    public static void setUpClass() {
        injector = Guice.createInjector(new ClothoTestModule(), new JongoModule());
        router = injector.getInstance(Router.class);
        SecurityManager securityManager = injector.getInstance(SecurityManager.class);
        SecurityUtils.setSecurityManager(securityManager);
        ClothoRealm realm = injector.getInstance(ClothoRealm.class);
        TestUtils.setupTestUsers(realm);
    }

    @AfterClass
    public static void tearDownClass() {
    }

    @Before
    public void setUp() {
        injector.getInstance(Persistor.class).deleteAll();
        ids = TestUtils.setupTestData(injector.getInstance(Persistor.class));
    }

    @After
    public void tearDown() {
        injector.getInstance(Persistor.class).deleteAll();
        ids = TestUtils.setupTestData(injector.getInstance(Persistor.class));
        
    }
    // TODO add test methods here.
    // The methods must be annotated with annotation @Test. For example:
    //
    // @Test
    // public void hello() {}

    @Test
    public void getAll() {
        TestConnection connection = new TestConnection("getTest");
        Message message = new Message(Channel.getAll, new String[]{"Test Part 1"} , "1");
        sendMessage(message, connection);
        Message returnMessage = connection.messages.get(1);
        assertMatch(message, returnMessage);
        assertEquals("Test Part 1", ((Map) ((List)returnMessage.data).get(0)).get("name").toString());
    }

    private void assertMatch(Message m1, Message m2) {
        assertEquals(m1.channel, m2.channel);
        assertEquals(m1.requestId, m2.requestId);
    }

    @Test
    public void createAll() {
        ObjectId id = new ObjectId();
        //XXX: this is actually bad data
        Map<String, Object> newPart = new HashMap<>();
        newPart.put("name", "Test Part 4");
        newPart.put("description", "previously unsaved test part");
        newPart.put("sequence", "CCCCCCCCCCCCCCCCCCCCCCCC");
        newPart.put("format", "FreeForm");
        newPart.put("id", id.toString());

        TestConnection connection = new TestConnection("createTest");
        Message message = new Message(Channel.createAll, new Map[]{newPart}, "2");
        sendMessage(message, connection);
        Message returnMessage = connection.messages.get(2);
        assertMatch(message, returnMessage);
        assertEquals(id.toString(), ((List)returnMessage.data).get(0).toString());
    }

    @Test
    public void query() {
        Map<String, Object> query = new HashMap<>();
        query.put("schema", "org.clothocad.model.Part");
        TestConnection connection = new TestConnection("queryTest");
        Message message = new Message(Channel.query, query, "3");
        sendMessage(message, connection);
        Message returnMessage = connection.messages.get(1);
        assertMatch(message, returnMessage);
        int size = ((List) returnMessage.data).size();
        assertNotEquals(0, size);
        //assertEquals(3, ((List) returnMessage.data).size());
        
        connection = new TestConnection("queryTest2");
        query = new HashMap<>();
        query.put("schema", "org.clothocad.model.BasicPart");
        message = new Message(Channel.query, query, "4");
        sendMessage(message, connection);
        returnMessage = connection.messages.get(1);
        assertMatch(message, returnMessage);
        assertNotEquals(0, ((List) returnMessage.data).size());
    }

    
    @Test
    public void destroyAll(){
        TestConnection connection = new TestConnection("destroyTest");
        List<String> stringIds = new ArrayList<>();
        for (ObjectId id : ids){
            stringIds.add(id.toString());
        }
        
        Message message = new Message(Channel.destroyAll, stringIds, "4");
        sendMessage(message, connection);
        Persistor persistor = injector.getInstance(Persistor.class);
        for (ObjectId id : ids){
            try {
                persistor.getAsJSON(id);
                fail();
            } catch (EntityNotFoundException e) {}
        }
    }
    
    @Test
    public void setAll(){
        TestConnection connection = new TestConnection("setTest");
        List<Map<String,Object>> specs = new ArrayList<>();
        
        for (ObjectId id : ids){
            Map<String, Object> spec = new HashMap<>();
            spec.put("id", id.toString());
            spec.put("name", "set");
            specs.add(spec);
        }
        
        Persistor persistor = injector.getInstance(Persistor.class);
        Message message = new Message(Channel.setAll, specs, "5");
        sendMessage(message, connection);
        for (ObjectId id : ids){
            Map<String,Object> result = persistor.getAsJSON(id);
            assertEquals("set",result.get("name"));
            //simple way to check that set has not disrupted the rest of the data in the object
            assertNotNull(result.get("schema"));
        }       
    }
    
    @Test
    public void testConstructFunction() {
        String script = 
                  "var data = {};\n"
                + "data.name = \"reverse1\";\n"
                + "data.language = \"JAVASCRIPT\";\n"
                + "data.schema = \"org.clothocad.core.datums.Function\";\n"
                + "data.code = \"function(sequence) { return sequence.split('').reverse().join('');};\";\n"
                + "data.arguments = [{name:'sequence', type:'String'}];\n"
                + "\n"
                + "clotho.create(data);\n"
                + "\n"
                + "clotho.run(\"reverse1\", [\"AAACCC\"]);";
        TestConnection connection = new TestConnection("constructFunction");
        
        sendMessage(new Message(Channel.submit, script, "6"), connection);
        
        assertEquals("CCCAAA", connection.messages.get(1).data);
        
    }

        
    //Test for persisting values in the scripting environment
    
    @Test
    public void mindPersistenceTest() {
        TestConnection connection = new TestConnection("persistenceTest");
        Map<String,String> credentials = new HashMap<>();
        credentials.put("username", "testuser");
        credentials.put("password", "password");
        
        //login as testuser
        sendMessage(new Message(Channel.login, credentials, "7"), connection);
        Object data = connection.messageDataByChannelAndId.get(Channel.login.name()+"7");
        assertTrue((Boolean) data);
        //set value
        sendMessage(new Message(Channel.submit, "var persistMe = 42 ", "8"), connection);
        //logout
        sendMessage(new Message(Channel.logout, "", "9"), connection);
        //check that value is not available to anonymous user
        sendMessage(new Message(Channel.submit, "persistMe", "9"), connection);
        data = connection.messageDataByChannelAndId.get(Channel.submit.name()+"9");
        assertNotEquals(data, 42);
        //login again as testuser
        sendMessage(new Message(Channel.login, credentials, "10"), connection);
        //check value is available again
        sendMessage(new Message(Channel.submit, "persistMe", "11"), connection);
        data = connection.messageDataByChannelAndId.get(Channel.submit.name()+"11");
        assertEquals(data, 42);
    }
    
    //TODO:
    public void testConstructFunctionStreamlined() {
        String script = 
                  "var reverse = function(sequence) {\n"
                + "	return sequence.split('').reverse().join('');\n"
                + "};\n"
                + "\n"
                + "clotho.save(reverse);\n"
                + "\n"
                + "reverse(\"AAACCC\");";
    }
    
    
    
    
    private void sendMessage(Message message, ClientConnection connection) {
        String stringMessage = JSON.serializeForExternal(message);
        try {
            message = JSON.mapper.readValue(stringMessage, Message.class);
        } catch (IOException ex) {
        }
        router.receiveMessage(connection, message);
    }
}
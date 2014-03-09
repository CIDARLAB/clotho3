package org.clothocad.core.communication;

import com.google.inject.Guice;
import com.google.inject.Injector;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.persistence.EntityNotFoundException;
import org.apache.shiro.mgt.DefaultSecurityManager;
import org.apache.shiro.SecurityUtils;
import org.bson.types.ObjectId;
import org.clothocad.core.communication.Channel;
import org.clothocad.core.communication.ClientConnection;
import org.clothocad.core.communication.Message;
import org.clothocad.core.communication.Router;
import org.clothocad.core.persistence.mongodb.MongoDBModule;
import org.clothocad.core.persistence.Persistor;
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
        injector = Guice.createInjector(new ClothoTestModule(), new MongoDBModule());
        router = injector.getInstance(Router.class);
        SecurityUtils.setSecurityManager(new DefaultSecurityManager());
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
    public void getAll() throws IOException {
        TestConnection connection = new TestConnection("getTest");
        final Message message = new Message(
            Channel.getAll,
            new String[] {"Test Part 1"},
            "1",
            null
        );
        sendMessage(message, connection);
        Message returnMessage = connection.messages.get(1);
        assertMatch(message, returnMessage);
        assertEquals("Test Part 1", ((Map) ((List)returnMessage.getData()).get(0)).get("name").toString());
    }

    private void assertMatch(Message m1, Message m2) {
        assertEquals(m1.getChannel(), m2.getChannel());
        assertEquals(m1.getRequestId(), m2.getRequestId());
    }

    @Test
    public void createAll() throws IOException {
        ObjectId id = new ObjectId();
        //XXX: this is actually bad data
        Map<String, Object> newPart = new HashMap<>();
        newPart.put("name", "Test Part 4");
        newPart.put("description", "previously unsaved test part");
        newPart.put("sequence", "CCCCCCCCCCCCCCCCCCCCCCCC");
        newPart.put("format", "FreeForm");
        newPart.put("id", id.toString());

        TestConnection connection = new TestConnection("createTest");
        final Message message = new Message(
            Channel.createAll,
            new Map[] {newPart},
            "2",
            null
        );
        sendMessage(message, connection);
        Message returnMessage = connection.messages.get(2);
        assertMatch(message, returnMessage);
        assertEquals(id.toString(), ((List)returnMessage.getData()).get(0).toString());
    }

    @Test
    public void query() throws IOException {
        Map<String, Object> query = new HashMap<>();
        query.put("schema", "Part");
        TestConnection connection = new TestConnection("queryTest");
        Message message = new Message(Channel.query, query, "3", null);
        sendMessage(message, connection);
        Message returnMessage = connection.messages.get(1);
        assertMatch(message, returnMessage);
        assertEquals(55, ((List) returnMessage.getData()).size());
        //assertEquals(3, ((List) returnMessage.data).size());
        
        connection = new TestConnection("queryTest2");
        query = new HashMap<>();
        query.put("schema", "BasicPart");
        message = new Message(Channel.query, query, "4", null);
        sendMessage(message, connection);
        returnMessage = connection.messages.get(1);
        assertMatch(message, returnMessage);
        assertEquals(54, ((List) returnMessage.getData()).size());
    }

    
    @Test
    public void destroyAll() throws IOException {
        TestConnection connection = new TestConnection("destroyTest");
        List<String> stringIds = new ArrayList<>();
        for (ObjectId id : ids){
            stringIds.add(id.toString());
        }
        
        final Message message = new Message(
            Channel.destroyAll,
            stringIds,
            "4",
            null
        );
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
    public void setAll() throws IOException {
        TestConnection connection = new TestConnection("setTest");
        List<Map<String,Object>> specs = new ArrayList<>();
        
        for (ObjectId id : ids){
            Map<String, Object> spec = new HashMap<>();
            spec.put("id", id.toString());
            spec.put("name", "set");
            specs.add(spec);
        }
        
        Persistor persistor = injector.getInstance(Persistor.class);
        sendMessage(
            new Message(
                Channel.setAll,
                specs,
                "5",
                null
            ),
            connection
        );
        for (ObjectId id : ids){
            Map<String,Object> result = persistor.getAsJSON(id);
            assertEquals("set",result.get("name"));
            assertNotNull(result.get("schema"));
        }       
    }
    
    @Test
    public void testConstructFunction() throws IOException {
        String script = 
                  "var data = {};\n"
                + "data.name = \"reverse1\";\n"
                + "data.language = \"JAVASCRIPT\";\n"
                + "data.schema = \"Function\";\n"
                + "data.code = \"function(sequence) { return sequence.split('').reverse().join('');};\";\n"
                + "data.arguments = [{name:'sequence', type:'String'}];\n"
                + "\n"
                + "clotho.create(data);\n"
                + "\n"
                + "clotho.run(\"reverse1\", [\"AAACCC\"]);";
        TestConnection connection = new TestConnection("constructFunction");
        
        sendMessage(
            new Message(
                Channel.submit,
                script,
                "6",
                null
            ),
            connection
        );
        
        assertEquals("CCCAAA", connection.messages.get(1).getData());
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
    
    
    
    
    private void sendMessage(Message message, ClientConnection connection)
    throws IOException {
        router.receiveMessage(
            connection,
            JSON.mapper.readValue(JSON.serialize(message), Message.class)
        );
    }
}

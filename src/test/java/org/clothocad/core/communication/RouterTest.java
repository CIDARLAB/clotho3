package org.clothocad.core.communication;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.persistence.EntityNotFoundException;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.util.AuthorizedShiroTest;
import org.clothocad.core.util.JSON;
import org.clothocad.core.util.TestUtils;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author spaige
 */
public class RouterTest extends AuthorizedShiroTest{

    private static Router router;
    private static List<ObjectId> ids;

    public RouterTest() {
            router = injector.getInstance(Router.class);

    }

    @BeforeClass
    public static void setUpClass() {

    }


    @Before
    public void setUp() {
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
            new String[] {ids.get(0).toString()},
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
        List ids = (List) connection.messageDataByChannelAndId.get(Channel.createAll.toString()+"2");
        assertEquals(id.toString(), ids.get(0).toString());
    }

    @Test
    public void query() throws IOException {
        Map<String, Object> query = new HashMap<>();
        query.put("schema", "org.clothocad.model.Part");
        TestConnection connection = new TestConnection("queryTest");
        Message message = new Message(Channel.query, query, "3", null);
        sendMessage(message, connection);
        Message returnMessage = connection.messages.get(1);
        assertMatch(message, returnMessage);
        int size = ((List) returnMessage.getData()).size();
        assertNotEquals(0, size);
        //assertEquals(3, ((List) returnMessage.data).size());
        
        connection = new TestConnection("queryTest2");
        query = new HashMap<>();
        query.put("schema", "org.clothocad.model.BasicPart");
        message = new Message(Channel.query, query, "4", null);
        sendMessage(message, connection);
        returnMessage = connection.messages.get(1);
        assertMatch(message, returnMessage);
        assertNotEquals(0, ((List) returnMessage.getData()).size());
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
                + "data.schema = \"org.clothocad.core.datums.Function\";\n"
                + "data.code = \"function(sequence) { return sequence.split('').reverse().join('');};\";\n"
                + "data.arguments = [{name:'sequence', type:'String'}];\n"
                + "\n"
                + "reverse1 = clotho.create(data);\n"
                + "\n"
                + "clotho.run(reverse1, [\"AAACCC\"]);";
        TestConnection connection = new TestConnection("constructFunction");
        
        Map submission = new HashMap();
        submission.put("query", script);
        submission.put("tokens", new ArrayList());

        sendMessage(new Message(Channel.submit, submission, "6"), connection);
        
        assertEquals("CCCAAA", connection.messages.get(1).getData());
        
    }

        
    //Test for persisting values in the scripting environment
    
    @Test
    public void mindPersistenceTest()  throws IOException {
        TestConnection connection = new TestConnection("persistenceTest");
        Map<String,String> credentials = new HashMap<>();
        credentials.put("username", "testuser");
        credentials.put("password", "password");
        
        Map submission = new HashMap();
        submission.put("tokens", new ArrayList());
        
        //login as testuser
        sendMessage(new Message(Channel.login, credentials, "7"), connection);
        Object data = connection.messageDataByChannelAndId.get(Channel.login.name()+"7");
        assertTrue((Boolean) data);
        //set value
        submission.put("query", "var persistMe = 42 ");
        sendMessage(new Message(Channel.submit, submission, "8"), connection);
        //logout
        sendMessage(new Message(Channel.logout, "", "9"), connection);
        //check that value is not available to anonymous user
        submission.put("query", "persistMe");
        sendMessage(new Message(Channel.submit, submission, "9"), connection);
        data = connection.messageDataByChannelAndId.get(Channel.submit.name()+"9");
        assertNotEquals(data, 42);
        //login again as testuser
        sendMessage(new Message(Channel.login, credentials, "10"), connection);
        //check value is available again
        sendMessage(new Message(Channel.submit, submission, "11"), connection);
        data = connection.messageDataByChannelAndId.get(Channel.submit.name()+"11");
        assertEquals(42, data);
    }
    
       @Test
    public void crossConnectionMindPersistenceTest()  throws IOException {
        TestConnection connection = new TestConnection("crossConnectionPersistenceTest");
        Map<String,String> credentials = new HashMap<>();
        credentials.put("username", "testuser");
        credentials.put("password", "password");
        
        Map submission = new HashMap();
        submission.put("tokens", new ArrayList());
        
        //login as testuser
        sendMessage(new Message(Channel.login, credentials, "7"), connection);
        Object data = connection.messageDataByChannelAndId.get(Channel.login.name()+"7");
        assertTrue((Boolean) data);
        //set value
        submission.put("query", "var persistMe = 42 ");
        sendMessage(new Message(Channel.submit, submission, "8"), connection);
        //logout
        sendMessage(new Message(Channel.logout, "", "9"), connection);
        //login again as testuser on different connection
        connection = new TestConnection("differentConnection");
        sendMessage(new Message(Channel.login, credentials, "10"), connection);
        //check value is available again
        sendMessage(new Message(Channel.submit, submission, "11"), connection);
        data = connection.messageDataByChannelAndId.get(Channel.submit.name()+"11");
        assertEquals(42, data);
    } 
    
    @Test
    public void scriptedGet() throws IOException{
        Map submission = new HashMap();
        submission.put("tokens", new ArrayList());
        
        TestConnection connection = new TestConnection("test");
        submission.put("query", "clotho.get(\"org.clothocad.trails.WritingTrails\")");
        sendMessage(new Message(Channel.submit, submission, "1"), connection);
        //make sure that sub-objects contain expected properties
        Map<String,Object> data = (Map) connection.messageDataByChannelAndId.get(Channel.submit.name()+"1");
        assertNotNull(data);
        assertTrue(data.containsKey("contents"));
        assertTrue(data.get("contents") instanceof List);
        assertTrue( ((Map) ((List) data.get("contents")).get(0)).containsKey("pages"));
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

    private void sendMessage(Message message, ClientConnection connection) throws IOException {
        String stringMessage = JSON.serializeForExternal(message);
        message = JSON.mapper.readValue(stringMessage, Message.class);
        router.receiveMessage(connection, message);
    }
}
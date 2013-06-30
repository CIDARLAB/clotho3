/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.testers.communication;

import com.fasterxml.jackson.core.JsonParseException;
import com.fasterxml.jackson.databind.JsonMappingException;
import com.google.inject.Guice;
import com.google.inject.Injector;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.bson.types.ObjectId;
import org.clothocad.core.ClothoModule;
import org.clothocad.core.layers.communication.Channel;
import org.clothocad.core.layers.communication.Message;
import org.clothocad.core.layers.communication.Router;
import org.clothocad.core.layers.communication.connection.ClientConnection;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.testers.MongoDBTestModule;
import org.clothocad.core.util.JSON;
import org.clothocad.core.utils.TestUtils;
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
        injector = Guice.createInjector(new ClothoModule(), new MongoDBTestModule());
        router = Router.get();
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
    }
    // TODO add test methods here.
    // The methods must be annotated with annotation @Test. For example:
    //
    // @Test
    // public void hello() {}

    @Test
    public void get() {
        TestConnection connection = new TestConnection("getTest");
        Message message = new Message(Channel.get, "Test Part 1", "1");
        sendMessage(message, connection);
        Message returnMessage = connection.messages.get(1);
        assertMatch(message, returnMessage);
        assertEquals("Test Part 1", ((Map) returnMessage.data).get("name").toString());
    }

    private void assertMatch(Message m1, Message m2) {
        assertEquals(m1.channel, m2.channel);
        assertEquals(m1.requestId, m2.requestId);
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

        TestConnection connection = new TestConnection("createTest");
        Message message = new Message(Channel.create, newPart, "2");
        sendMessage(message, connection);
        Message returnMessage = connection.messages.get(1);
        assertMatch(message, returnMessage);
        assertEquals(id.toString(), returnMessage.data.toString());
    }

    @Test
    public void query() {
        Map<String, Object> query = new HashMap<>();
        query.put("schema", "Part");
        TestConnection connection = new TestConnection("queryTest");
        Message message = new Message(Channel.query, query, "3");
        sendMessage(message, connection);
        Message returnMessage = connection.messages.get(1);
        assertMatch(message, returnMessage);
        assertEquals(3, ((List) returnMessage.data).size());
    }

    private void sendMessage(Message message, ClientConnection connection) {
        String stringMessage = JSON.serialize(message);
        try {
            message = JSON.mapper.readValue(stringMessage, Message.class);
        } catch (IOException ex) {
        }
        router.receiveMessage(connection, message);
    }
}
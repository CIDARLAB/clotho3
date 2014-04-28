/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.communication;

import com.google.inject.Guice;
import com.google.inject.Injector;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.apache.shiro.SecurityUtils;
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
import static org.junit.Assert.*;
import org.junit.Test;
import org.python.google.common.collect.Lists;

/**
 *
 * @author spaige
 */
public class MessageOptionsTest {
    private static Router router;
    private static Injector injector;
    private static List<ObjectId> ids;
    
    @BeforeClass
    public static void setUpClass() {
        injector = Guice.createInjector(new ClothoTestModule(), new JongoModule());
        router = injector.getInstance(Router.class);
        org.apache.shiro.mgt.SecurityManager securityManager = injector.getInstance(org.apache.shiro.mgt.SecurityManager.class);
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
    
    private void sendMessage(Message message, ClientConnection connection) throws IOException {
        String stringMessage = JSON.serializeForExternal(message);
        message = JSON.mapper.readValue(stringMessage, Message.class);
        router.receiveMessage(connection, message);
    }
    
    @Test
    public void testMaxResults() throws IOException{
        TestConnection conn = new TestConnection("maxResults");
        Map<String,Object> data = new HashMap<>();
        data.put("schema", "org.clothocad.model.Part");
        Map<MessageOption,Object> options = new HashMap<>();
        options.put(MessageOption.maxResults, 1);
        final Message message = new Message(Channel.query, data, "1", options);
        sendMessage(message,conn);
        
        List results = (List) conn.messageDataByChannelAndId.get(Channel.query.name()+"1");
        assertEquals(1,results.size());
    }
    
    @Test
    public void testMute() throws IOException{
        TestConnection conn = new TestConnection("mute");
        Map<String,Object> data = new HashMap<>();
        data.put("schema", "org.clothocad.model.Part");
        Map<MessageOption,Object> options = new HashMap<>();
        options.put(MessageOption.mute, true);
        final Message message = new Message(Channel.query, data, "1", options);
        sendMessage(message,conn);
        
        for (Message clientMessage : conn.messages){
            assertNotEquals(Channel.say,clientMessage.getChannel());
        }
    }
    
    @Test
    public void testFilter() throws IOException{
        TestConnection conn = new TestConnection("mute");
        Map<String,Object> data = new HashMap<>();
        data.put("schema", "org.clothocad.model.Part");
        Map<MessageOption,Object> options = new HashMap<>();
        options.put(MessageOption.filter, Lists.newArrayList("name"));
        final Message message = new Message(Channel.query, data, "1", options);
        sendMessage(message,conn);
        
        List<Map<String,Object>> results = (List) conn.messageDataByChannelAndId.get(Channel.query.name()+"1");
        for (Map<String,Object> result : results){
            assertEquals(2, result.keySet().size());
            assertTrue(result.containsKey("name"));
            assertTrue(result.containsKey("id"));
        }
    }
}

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.communication;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import static org.junit.Assert.*;
import org.junit.Test;
import org.python.google.common.collect.Lists;

/**
 *
 * @author spaige
 */
public class MessageOptionsTest extends AbstractRouterTest{
    public MessageOptionsTest() {
        super();
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
    public void testMuteRunModuleMethod() throws IOException{
        TestConnection conn = new TestConnection("mute");
        Map<String,Object> data = new HashMap<>();
        data.put("id", "org.clothocad.test.testModule2");
        data.put("function", "moduleMethod");
        data.put("args", new String[]{});
        Map<MessageOption,Object> options = new HashMap<>();
        options.put(MessageOption.mute, true);
        final Message message = new Message(Channel.run , data, "1", options);
        sendMessage(message,conn);
        
        for (Message clientMessage : conn.messages){
            assertNotEquals(Channel.say,clientMessage.getChannel());
        }
        assertEquals(2, conn.messageDataByChannelAndId.get("run1"));
    }
    
    @Test
    public void testMuteRunFunctionMethod() throws IOException{
        TestConnection conn = new TestConnection("mute");
        Map<String,Object> data = new HashMap<>();
        data.put("id", "org.clothocad.test.moduleTestFunction");
        data.put("args", new Integer[]{1});
        Map<MessageOption,Object> options = new HashMap<>();
        options.put(MessageOption.mute, true);
        final Message message = new Message(Channel.run , data, "1", options);
        sendMessage(message,conn);
        
        for (Message clientMessage : conn.messages){
            assertNotEquals(Channel.say,clientMessage.getChannel());
        }
        assertEquals(4.0, conn.messageDataByChannelAndId.get("run1"));
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

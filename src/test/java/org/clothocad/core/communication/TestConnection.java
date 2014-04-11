/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.communication;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author spaige
 */
public class TestConnection extends ClientConnection{

    public TestConnection(String id) {
        super(id);
    }

    public List<Message> messages = new ArrayList<>();
    public Map<String, Object> messageDataByChannelAndId= new HashMap<>();
    
    @Override
    public void send(Message msg) {
        messages.add(msg);
        messageDataByChannelAndId.put(msg.getChannel().name() + msg.getRequestId(), msg.getData());
    }
    
}

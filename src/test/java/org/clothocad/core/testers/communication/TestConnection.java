/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.testers.communication;

import java.util.ArrayList;
import java.util.List;
import org.clothocad.core.layers.communication.Message;
import org.clothocad.core.layers.communication.connection.ClientConnection;

/**
 *
 * @author spaige
 */
public class TestConnection extends ClientConnection{

    public TestConnection(String id) {
        super(id);
    }

    public List<Message> messages = new ArrayList<>();
    
    @Override
    public void send(Message msg) {
        messages.add(msg);
    }
    
}

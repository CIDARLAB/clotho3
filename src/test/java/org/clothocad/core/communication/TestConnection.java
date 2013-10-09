/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.communication;

import java.util.ArrayList;
import java.util.List;
import org.clothocad.core.communication.Message;
import org.clothocad.core.communication.ClientConnection;

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

package org.clothocad.core.communication;

import org.clothocad.core.communication.Channel;
import org.clothocad.core.communication.Message;

public abstract class ClientConnection {

    private String id;

    public ClientConnection(String id) {
        this.id = id;
    }

    public String getId() {
        return this.id;
    }

    public abstract void send(Message msg);

    //TODO (make abstract)
    public void deregister(Channel channel, String requestId) {
    }
}

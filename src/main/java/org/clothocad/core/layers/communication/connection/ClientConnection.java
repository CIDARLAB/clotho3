package org.clothocad.core.layers.communication.connection;

import org.clothocad.core.layers.communication.Channel;
import org.clothocad.core.layers.communication.Message;

public abstract class ClientConnection {

	private String id;
	public ClientConnection(String id) {
		this.id = id;
	}
	public String getId() {
		return this.id;
	}
        
        public abstract void send(Message msg);

    public void deregister(Channel channel, String requestId) {
    }
}

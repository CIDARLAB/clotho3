package org.clothocad.core.layers.communication.connection.ws;

import java.io.IOException;

import org.clothocad.core.layers.communication.ClothoConstants;
import org.clothocad.core.layers.communication.Router;
import org.clothocad.core.layers.communication.connection.ClientConnection;
import org.eclipse.jetty.websocket.WebSocket;
import org.json.JSONObject;

public class ClothoWebSocket 
	extends ClientConnection 
	implements WebSocket.OnTextMessage {

	private WebSocket.Connection connection;
	
	public ClothoWebSocket(String id) {
            super(id);
	}
	
	@Override
	public void onClose(int closeCode, String message) {
	}

	public void sendMessage(String data) 
                throws IOException {
            if(connection.isOpen()) {
                connection.sendMessage(data);
            }
	}

	@Override
	public void onMessage(String message) {
		// here, we need to forward the message dependent on its content

            JSONObject json = null;
            try {
                json = new JSONObject(message);
            } catch(Exception e) {
                // invalid json object -> do some error handling
                e.printStackTrace();
            }					
		
            // do the unmarshaling
            if(null != json) {
                try {					
                    // do some routing
                    Router.get().receiveMessage(
                        this, 
                        json.getString(ClothoConstants.CHANNEL), 
                        json);
                } catch(Exception e) {
                    e.printStackTrace();
                }
            }
	}

	public boolean isOpen() {
		return connection.isOpen();
	}

	@Override
	public void onOpen(Connection connection) {
		this.connection = connection;
                //Close out after 1 hour idle time
                //connection.setMaxIdleTime(3600000); 

		//WebSocketTable.put(this.getId(), this);
		
		// TODO: store the connection information into the Mind
		
	}
}

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

	private Connection connection;
	
	public ClothoWebSocket(String id) {
		super(id);
	}
	
	@Override
	public void onClose(int closeCode, String message) {
	}

	public void sendMessage(String data) throws IOException {
		System.out.println("[ClothoWebSocket.sendMessage] -> "+data);		
		connection.sendMessage(data);
	}

	@Override
	public void onMessage(String message) {
		// here, we need to forward the message dependent on its content
		
		System.out.println("[ClothoWebSocket.onMessage] -> "+message);

		/**
		obj = JSON.parse(obj);
		var channel = obj.channel;
		var data = obj.data;
		
		//note - channel reserved for serverAPI
		if (channel == "$clotho") {
		    console.log("SOCKET\tchannel $clotho");
		    internal_trigger(data.channel, data.data);
		    return;
		}
		
		// it's the ClientAPI method's responsibility to handle data appropriately.
		if (typeof ClientAPI[channel] == 'function') {
		    console.log("SOCKET\tmapping to ClientAPI - " + channel);
		    ClientAPI[channel](data);
		}
		//for custom listeners attached
		else if (typeof customChannels[channel] == 'function') {
		    console.log("SOCKET\tmapping to custom listeners - " + channel);
		    trigger(channel, data);
		}
		// don't know what to do, so publish to PubSub
		else {
		    console.log("SOCKET\tno listener found for channel: " + channel);
		    PubSub.publish(channel, data);
		}
		 **/
		
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

		//WebSocketTable.put(this.getId(), this);
		
		// TODO: store the connection information into the Mind
		
	}
}

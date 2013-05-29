package org.clothocad.core.layers.communication;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.clothocad.core.layers.communication.activemq.ClothoMessageProducer;
import org.clothocad.core.layers.communication.connection.ClientConnection;
import org.clothocad.core.layers.communication.connection.apollo.ApolloConnection;
import org.clothocad.core.layers.communication.connection.ws.ClothoWebSocket;
import org.clothocad.core.layers.communication.protocol.ActionType;
import org.json.JSONException;
import org.json.JSONObject;

public class Router {

	/** DOUBLE-CHECKED LOCKING **/
	private static volatile Router router;

    public static Router get() {
    	Router result = router;
		if(result == null) {
			synchronized(Router.class) {
				result = router;
				if(result == null) {
					router = result = new Router();
				}
			}
		}
		return result;
    }
    
    private ExecutorService pool;
    
    public Router() {    	
    	
    	// create the thread-pool of ssAPI objects
    	this.pool = Executors.newFixedThreadPool(100);
    }
    
    // send message    
	public void sendMessage(
			ClientConnection connection, String channel, JSONObject data) {
		
		
		// create the message's JSON object
		JSONObject message = new JSONObject();
		try {
			message.put(ClothoConstants.CHANNEL, channel);
			message.put(ClothoConstants.DATA, data);

			if(connection instanceof ClothoWebSocket) {
				ClothoWebSocket websocket = (ClothoWebSocket)connection;
				if(null != websocket) {				
					websocket.sendMessage(message.toString());
				}
			} else if(connection instanceof ApolloConnection) {				
				new ClothoMessageProducer((ApolloConnection)connection).onSuccess(message);
			}
			
		} catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	// receive message
	public void receiveMessage(ClientConnection connection, String channel, JSONObject json) {
		System.out.println("[Router.receiveMessage] -> "+connection+", "+channel+", "+json.toString());
		
		try {
			if(Channel.ACCESS.toString().equalsIgnoreCase(channel)) {
				// we immediately respond to the client on the synchronous ACCESS channel
				
				JSONObject response = this.access(json);
				
				this.sendMessage(
						connection, channel, response);
			} else if(Channel.EXECUTION.toString().equalsIgnoreCase(channel)) {								
				
				// here we should get rid of giving socket_id, channel, and the JSON object to the Executor
				this.execute(connection, json);
				
				// ultimately:
				// this.execute(json);
				
			} else if(Channel.NOTIFICATION.toString().equalsIgnoreCase(channel)) {
				// TODO
			} else if(Channel.AUTOCOMPLETION.toString().equalsIgnoreCase(channel)) {
				// TODO

			} else if("say".equals(channel)) {
			}
		} catch(Exception e) {
			e.printStackTrace();
		}
	}	
	
	private JSONObject access(JSONObject json) 
			throws JSONException {

		// here, we need to invoke the data-base layer to access to requested data
		JSONObject response = new JSONObject();
		
		// first, get the action of the JSON object
		String sAction = json.getString("action");

		if(ActionType.LOGIN.toString().equals(sAction)) {			

			// here, we need to call the database layer to check if the login-data is valid
			/**
			boolean b = Database.get().login(
					String.valueOf(datum.get("username")), 
					String.valueOf(datum.get("password")));
			**/
			
			response.put("valid", true);
		}

		return response;
	}
	
	private void execute(ClientConnection connection, JSONObject json) 
			throws Exception {
		// first, get the action of the JSON object
		String action = json.getString(ClothoConstants.ACTION);

		// get the correlation-id too
		String correlationId = json.getString(ClothoConstants.CORRELATION_ID);
		
		// and the authentication key
		String authKey = json.getString(ClothoConstants.AUTHENTICATION);
		
		// that's the data of the incoming request
		
		JSONObject data = new JSONObject(json.get(ClothoConstants.DATA));
		
		this.pool.execute(
				new RequestProcessor(connection, action, correlationId, data));		
	}
}

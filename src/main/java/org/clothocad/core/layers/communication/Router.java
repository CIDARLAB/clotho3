package org.clothocad.core.layers.communication;

import org.clothocad.core.aspects.Persistor;
import org.clothocad.core.layers.communication.activemq.ClothoMessageProducer;
import org.clothocad.core.layers.communication.activemq.ClothoPublisher;
import org.clothocad.core.layers.communication.protocol.ActionType;
import org.clothocad.core.layers.execution.Executor;
import org.clothocad.core.layers.execution.ExecutorThread;
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
    
	// send message
	public void sendMessage(String socket_id, String channel, JSONObject json) {
		try {
			new ClothoMessageProducer(socket_id).onSuccess(json);
		} catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	// publish messages to subscribed clients
	public void publish(JSONObject json) {
		try {
			new ClothoPublisher().publish(json);
		} catch(Exception e) {}
	}
	
	// receive message
	public void receiveMessage(String socket_id, String channel, JSONObject json) {
		try {
			if(Channel.ACCESS.toString().equals(channel)) {
				// we immediately respond to the client on the synchronous ACCESS channel
				this.sendMessage(
						socket_id,  
						channel, 
						this.access(json));
			} else if(Channel.EXECUTION.toString().equals(channel)) {
				// here we should get rid of giving socket_id, channel, and the JSON object to the Executor
				//this.execute(json);
				new Thread(new ExecutorThread(socket_id, channel, json)).start();				
			} else if(Channel.NOTIFICATION.toString().equals(channel)) {
				// TODO
			} else if(Channel.AUTOCOMPLETION.toString().equals(channel)) {
				// TODO
			}
		} catch(Exception e) {
			e.printStackTrace();
		}
	}	
	
	private JSONObject access(JSONObject json) 
			throws JSONException {

		// here, we need to invoke the data-base layer to access to requested data
		JSONObject ret = new JSONObject();
		
		// first, get the action of the JSON object
		String sAction = json.getString("action");

		if(ActionType.LOGIN.toString().equals(sAction)) {			

			// here, we need to call the database layer to check if the login-data is valid
			/**
			boolean b = Database.get().login(
					String.valueOf(datum.get("username")), 
					String.valueOf(datum.get("password")));
			**/
			
			ret.put("valid", true);
			
		} else if(ActionType.GET.toString().equals(sAction)) {
			/**
			Object obj = Database.get().get(
					String.valueOf(json.get("id")));
			**/
			
			//ret.put("datum", null);
		}

		return ret;
	}
}

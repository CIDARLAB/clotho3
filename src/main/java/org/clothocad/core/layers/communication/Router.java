package org.clothocad.core.layers.communication;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import org.bson.types.ObjectId;
import org.clothocad.core.aspects.Interpreter.AutoComplete;
import org.clothocad.core.aspects.Persistor;
import org.clothocad.core.datums.Doo;

import org.clothocad.core.layers.communication.activemq.ClothoMessageProducer;
import org.clothocad.core.layers.communication.connection.ClientConnection;
import org.clothocad.core.layers.communication.connection.apollo.ApolloConnection;
import org.clothocad.core.layers.communication.connection.ws.ClothoWebSocket;
import org.clothocad.core.layers.communication.mind.Mind;
import org.clothocad.core.layers.communication.protocol.ActionType;
import org.clothocad.model.Institution;
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
			ClientConnection connection, JSONObject message) {

		try {

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
//		System.out.println("[Router.receiveMessage] -> "+connection+", "+channel+", "+json.toString());
            
		try {
                    RouterDoo doo = new RouterDoo();
                    
                    String auth_key = null;
                    try {
                         auth_key = json.getString("auth_key");
                    } catch(Exception error) {
                    }
                    
                    Mind mind = Communicator.get().getMind(auth_key);
                    mind.setClientConnection(connection);
                    doo.mindId = mind.getUUID();
                    doo.message = json;
                    ServerSideAPI api = mind.getAPI();

			if(Channel.autocomplete.toString().equalsIgnoreCase(channel)) {
                            api.autocomplete(json.getJSONObject("data").getString("query"));
			} else if(Channel.submit.toString().equalsIgnoreCase(channel)) {								
                            api.submit(json.getJSONObject("data").getString("query"));
			} else if(Channel.login.toString().equalsIgnoreCase(channel)) {

			} else if(Channel.autocompleteDetail.toString().equalsIgnoreCase(channel)) {
                            String uuid = json.getString("data");
                            api.autocompleteDetail(uuid);
			} else if(Channel.logout.toString().equalsIgnoreCase(channel)) {

			} else if(Channel.changePassword.toString().equalsIgnoreCase(channel)) {								

			} else if(Channel.say.toString().equalsIgnoreCase(channel)) {
                            api.say(json.getJSONObject("query").getString("text"));
			} else if(Channel.log.toString().equalsIgnoreCase(channel)) {

			}
                        //Etcetera for remaining ssAPI methods
		} catch(Exception e) {
			e.printStackTrace();
		}
	}	
        
    /**
     * The Doo's that manage any Client-derived message.  Doo's handle even key
     * commands to avoid synchronization issues.
     */
    public class RouterDoo extends Doo {
        public RouterDoo() {
            super(null, false);
        }

        JSONObject message;
        ObjectId mindId;
    }
    
    
    
}

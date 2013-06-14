package org.clothocad.core.layers.communication;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import org.bson.types.ObjectId;
import org.clothocad.core.datums.Doo;
import static org.clothocad.core.layers.communication.Channel.autocompleteDetail;
import static org.clothocad.core.layers.communication.Channel.create;
import static org.clothocad.core.layers.communication.Channel.destroy;
import static org.clothocad.core.layers.communication.Channel.log;
import static org.clothocad.core.layers.communication.Channel.submit;

import org.clothocad.core.layers.communication.activemq.ClothoMessageProducer;
import org.clothocad.core.layers.communication.connection.ClientConnection;
import org.clothocad.core.layers.communication.connection.apollo.ApolloConnection;
import org.clothocad.core.layers.communication.connection.ws.ClothoWebSocket;
import org.clothocad.core.layers.communication.mind.Mind;
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
                    
                    Channel chanEnum = Channel.valueOf(channel);
                    String data = json.getString("data");
                    
                    switch(chanEnum) {
                        case autocomplete:
                            api.autocomplete(data);
                            break;
                        case autocompleteDetail:
                            api.autocompleteDetail(data);
                            break;
                        case submit:
                            api.submit(data);
                            break;
                        case clear:
                            api.clear();
                            break;
                        case login:
                            api.login(data, null);
                            break;
                        case logout:
                            api.logout();
                            break;
                        case changePassword:
                            api.changePassword(data);
                            break;
                        case learn:
                            api.learn(data, null);
                            break;

                        
                        case log:
                            api.log(data);
                            break;
                        case say:
                            api.say(data);
                            break;
                        case note:
                            api.note(data);
                            break;
                        case alert:
                            api.alert(data);
                            break;

                        case get:
                            api.get(data);
                            break;
                        case set:
                            api.set(data);
                            break;
                        case create:
                            api.create(data);
                            break;
                        case destroy:
                            api.destroy(data);
                            break;
                        case query:
                            api.query(data);
                            break;

                        case run:
                            api.run(data, null);
                            break;
                        case show:
                            api.show(data, null, null, null);
                            break;
                        case startTrail:
                            api.startTrail(data);
                            break;
                        case edit:
                            api.edit(data);
                            break;
                        case listen:
                            api.listen(data);
                            break;
                        case unlisten:
                            api.unlisten(data);
                            break;
                        default:
                            break;
                    }
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

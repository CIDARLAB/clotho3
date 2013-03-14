package org.clothocad.core.layers.communication.activemq;

import javax.jms.Message;
import javax.jms.Session;

import org.clothocad.core.layers.communication.Router;
import org.json.JSONObject;

public class MessageProcessor 
	implements Runnable {
	
	private Session session;
	private Message request;
	
	public MessageProcessor(Session session, Message request) {
		this.session = session;
		this.request = request;
	}
	
	public void run() {
		try {
			JSONObject json = new JSONObject(request.getStringProperty("request"));
	
			// create a new callback-handler and put it into 
			// a hashtable for later execution...
			String socket_id = CallbackHandlerTable.put(
					new CallbackHandler(this.session, request));
			
	    	// here, we utilize the Router to process the incoming JSON packet
	    	Router.get().receiveMessage(
	    			socket_id,
	    			json.getString("channel"), 
	    			json);       // TODO: improve this... 
		} catch(Exception e) {
			e.printStackTrace();
		}
	} 
	
	private static MessageProcessor msgProc = null;
	
	public static MessageProcessor get(Session session, Message request) {
		msgProc = new MessageProcessor(session, request);
		return msgProc;
	}
}

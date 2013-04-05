package org.clothocad.core.layers.communication.activemq;

import javax.jms.Message;
import javax.jms.Session;
import javax.jms.TextMessage;

import org.apache.activemq.transport.stomp.StompFrame;
import org.clothocad.core.layers.communication.Router;
import org.json.JSONObject;

public class MessageProcessor 
	implements Runnable {
	
	private Session session;
	private Message message;
	
	public MessageProcessor(Session session, Message message) {
		this.session = session;
		this.message = message;
	}
	
	public void run() {
		try {
			if(message.propertyExists("request")) {
				JSONObject json = new JSONObject(
						message.getStringProperty("request"));
		
				// create a new callback-handler and put it into 
				// a hashtable for later execution...
				String socket_id = CallbackHandlerTable.put(
						new CallbackHandler(this.session, message));
				
		    	// here, we utilize the Router to route the incoming JSON packet
		    	Router.get().receiveMessage(
		    			socket_id,
		    			json.getString("channel"), 
		    			json);       // TODO: improve this...
			}
			
			StompFrame frame = (StompFrame)message;
			System.out.println(frame.getBody());
			
			
			// INVALID REQUEST 
			// WHAT TO DO?
			//System.out.println(json);
		} catch(Exception e) {
			e.printStackTrace();
		}
	} 
	
	/**
	private static MessageProcessor msgProc = null;
	
	public static MessageProcessor get(Session session, Message request) {
		msgProc = new MessageProcessor(session, request);
		return msgProc;
	}
	**/
}

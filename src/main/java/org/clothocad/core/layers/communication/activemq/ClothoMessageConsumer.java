package org.clothocad.core.layers.communication.activemq;

import java.util.UUID;

import javax.jms.JMSException;
import javax.jms.Message;
import javax.jms.Session;

import org.clothocad.core.layers.communication.Router;
import org.fusesource.stomp.jms.message.StompJmsMessage;
import org.json.JSONException;
import org.json.JSONObject;

public class ClothoMessageConsumer 
	implements Runnable {
	
	private Session session;
	private StompJmsMessage message;
	
	public ClothoMessageConsumer(Session session, Message message) 
			throws Exception {
		this.session = session;
		
		if(message!=null && message instanceof StompJmsMessage) {			
			this.message = (StompJmsMessage)message;
		} else {
			throw new Exception("INVALID MESSAGE!");
		}
	}
	
	@Override
	public void run() {
		try {
			if(this.message.propertyExists("request")) {
				
				JSONObject json = new JSONObject(
						message.getStringProperty("request"));
				
				// get the message's correlation id
				String sCorrelationID = this.message.getJMSCorrelationID();
				if(null == sCorrelationID) {
					sCorrelationID = UUID.randomUUID().toString();
					this.message.setJMSCorrelationID(sCorrelationID);
				}
				
				// store the callback-handler in the callback-handler table
				CallbackHandlerTable.put(
						sCorrelationID, 
						new CallbackHandler(this.session, this.message));
				
				// route the message
				Router.get().receiveMessage(
						sCorrelationID,
						json.getString("channel"), 
						json);
			}
		} catch (JMSException e) {
			e.printStackTrace();
		} catch (JSONException e) {
			e.printStackTrace();
		}
	} 
}

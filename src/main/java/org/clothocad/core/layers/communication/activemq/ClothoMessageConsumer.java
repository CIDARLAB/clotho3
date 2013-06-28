package org.clothocad.core.layers.communication.activemq;

import com.fasterxml.jackson.core.JsonParseException;
import java.util.Map;
import java.util.UUID;

import javax.jms.JMSException;
import javax.jms.Session;
import org.clothocad.core.layers.communication.Channel;

import org.clothocad.core.layers.communication.Router;
import org.clothocad.core.layers.communication.connection.apollo.ApolloConnection;
import org.clothocad.core.util.JSON;
import org.fusesource.stomp.jms.message.StompJmsMessage;

public class ClothoMessageConsumer 
	implements Runnable {
	
	private Session session;
	private StompJmsMessage message;
	
	public ClothoMessageConsumer(Session session, javax.jms.Message message) 
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
				
				Map<String, Object> json = JSON.deserializeObject(message.getStringProperty("request"));
				
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
				
				ApolloConnection connection = new ApolloConnection(
						session, sCorrelationID);
				
				// route the message
				Router.get().receiveMessage(
						connection, new org.clothocad.core.layers.communication.Message(Channel.valueOf(json.get("channel").toString()),json));
			}
		} catch (JMSException e) {
			e.printStackTrace();
		} catch (JsonParseException e) {
			e.printStackTrace();
		}
	} 
}

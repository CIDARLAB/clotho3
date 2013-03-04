package org.clothocad.core.layers.communication.activemq;

import javax.jms.DeliveryMode;
import javax.jms.Message;
import javax.jms.MessageProducer;
import javax.jms.Session;

import org.clothocad.core.layers.communication.protocol.ActionType;
import org.json.JSONObject;

public class CallbackHandler {

	private Session session;
	private Message request;
	private MessageProducer replyProducer;
	
	public CallbackHandler(Session session, Message request) {
		this.session = session;
		this.request = request;
		
		try {
			// Setup a message producer to respond to messages from clients, we will get the destination
	        // to send to from the JMSReplyTo header field from a Message
	        this.replyProducer = this.session.createProducer(null);
	        this.replyProducer.setDeliveryMode(DeliveryMode.NON_PERSISTENT);
		} catch(Exception e) {
			e.printStackTrace();
		}
    }
	
	public void respond(JSONObject json) {
		try {
			Message response = session.createMessage();
			response.setJMSCorrelationID(request.getJMSCorrelationID());
			response.setStringProperty(ActionType.SERVER_RESPONSE.toString(), json.toString());
			
			replyProducer.send(
					request.getJMSReplyTo(), 
					response);			
		} catch(javax.jms.InvalidDestinationException e) {
			// something went wrong while responding to the client...
			System.err.println("The client does not run anymore...");
		} catch(javax.jms.JMSException e) {
			// something went wrong with the message ...
		}
	}
}

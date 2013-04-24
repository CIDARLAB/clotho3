package org.clothocad.core.layers.communication.activemq;

import javax.jms.DeliveryMode;
import javax.jms.Message;
import javax.jms.MessageProducer;
import javax.jms.Session;

import org.clothocad.core.layers.communication.protocol.ActionType;
import org.fusesource.stomp.jms.StompJmsDestination;
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
	        this.replyProducer = this.session.createProducer(
	        		new StompJmsDestination("/queue/CLOTHORESPONSE"));
	        this.replyProducer.setDeliveryMode(DeliveryMode.NON_PERSISTENT);

		} catch(Exception e) {
			e.printStackTrace();
		}
    }
	
	public void respond(JSONObject json) {
		try {
			Message response = session.createMessage();
			response.setJMSCorrelationID(request.getJMSCorrelationID());
			response.setStringProperty(ActionType.RESPONSE.toString(), json.toString());

			this.replyProducer.send(response);
		} catch(javax.jms.InvalidDestinationException e) {
			// something went wrong while responding to the client...
			e.printStackTrace();
		} catch(javax.jms.JMSException e) {
			// something went wrong with the message ...
		}
	}
}

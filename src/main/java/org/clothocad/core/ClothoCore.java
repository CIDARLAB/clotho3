package org.clothocad.core;

import javax.jms.Connection;
import javax.jms.JMSException;
import javax.jms.Message;
import javax.jms.MessageConsumer;
import javax.jms.MessageListener;
import javax.jms.MessageProducer;
import javax.jms.Session;

import org.clothocad.core.layers.communication.activemq.ClothoMessageConsumer;
import org.fusesource.stomp.jms.StompJmsConnectionFactory;
import org.fusesource.stomp.jms.StompJmsDestination;

public class ClothoCore 
		implements MessageListener {
	
	private Connection connection = null;
	private Session session = null;
	
	public ClothoCore() 
			throws Exception {
		StompJmsConnectionFactory factory = 
				new StompJmsConnectionFactory();
		factory.setBrokerURI("tcp://localhost:61613");
		
		this.connection = factory.createConnection("admin", "password");
		this.connection.start();
	}
	
	public void start() 
			throws Exception {
		
		// third, we create a session for the connection
        this.session = connection.createSession(
        		false, Session.CLIENT_ACKNOWLEDGE);

        // fourth, we create a consumer which reads incoming messages
        MessageConsumer consumer = session.createConsumer(
        		new StompJmsDestination("/queue/CLOTHO"));
        
        MessageProducer producer = session.createProducer(
        		new StompJmsDestination("/queue/CLOTHORESPONSE"));

        // finally, we listen for incoming messages
        consumer.setMessageListener(this);  
        		// for every incoming message, 
        		// the onMessage method gets invoked

        System.out.println("The Core is running...");
    }
	
    // the onMessage method gets executed if a message arrives on the given channel...
	@Override
	public void onMessage(Message message) {
		try {
			new Thread(
					new ClothoMessageConsumer(
							this.session, 
							message)).start();
		} catch (Exception e) {
			// currently, we ignore invalid messages
			e.printStackTrace();
		}
	}
	
	
	public void shutdown() 
			throws JMSException {
		if(null != connection) {
	        connection.close();
		}
	}
}

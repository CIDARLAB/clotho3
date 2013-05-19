package org.clothocad.ws;

import javax.jms.Connection;
import javax.jms.DeliveryMode;
import javax.jms.JMSException;
import javax.jms.Session;

import org.fusesource.stomp.jms.StompJmsConnectionFactory;
import org.fusesource.stomp.jms.StompJmsDestination;
import javax.jms.MessageProducer;


public class ApolloWSServer {

	public static void main(String[] args) {
        StompJmsConnectionFactory factory = new StompJmsConnectionFactory();
        factory.setBrokerURI("ws://localhost:61623");

        try {

            Connection connection = factory.createConnection("admin", "password");
            connection.start();
            
            Session session = connection.createSession(false, Session.AUTO_ACKNOWLEDGE);
            MessageProducer messageProducer = session.createProducer(
    				new StompJmsDestination("/queue/CLOTHO"));

    		messageProducer.setTimeToLive(10000);
    		messageProducer.setDeliveryMode(DeliveryMode.NON_PERSISTENT);
    		
            //responseConsumer = session.createConsumer(
            //		new StompJmsDestination("/queue/CLOTHORESPONSE"));
            
        } catch (JMSException e) {
        	e.printStackTrace();
        	// if there is an exception, it's highly possible that the server is not running
        	// -> go to ``offline'' mode
        }
	}

}

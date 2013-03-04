package org.clothocad.core.layers.communication.activemq;

import javax.jms.Connection;
import javax.jms.Destination;
import javax.jms.JMSException;
import javax.jms.Message;
import javax.jms.MessageConsumer;
import javax.jms.MessageListener;
import javax.jms.Session;

import org.apache.activemq.ActiveMQConnectionFactory;


/*
 * A channel listener listens for incoming messages on the channel/queue...
 * it's onMessage() method gets executed for every incoming message... 
 * and for incoming message, we create a thread which does the message processing... 
 */
public class ChannelListener 
		implements MessageListener {
	
    private int ackMode;
    private Session session;    
    private boolean transacted;

    public ChannelListener(String messageBrokerUrl) {
    	this.transacted = false;
        this.ackMode = Session.AUTO_ACKNOWLEDGE;
        
        ActiveMQConnectionFactory connectionFactory = 
        		new ActiveMQConnectionFactory(messageBrokerUrl);
        Connection connection;
        try {
            connection = connectionFactory.createConnection();
            connection.start();
            
            this.session = connection.createSession(this.transacted, ackMode);

            // this is the queue, the listener is listening to... 
            Destination adminQueue = this.session.createQueue("CLOTHO");
           
            // create a message consumer that listens on the queue
            MessageConsumer consumer = this.session.createConsumer(adminQueue);
            consumer.setMessageListener(this);
            
        } catch (JMSException e) {
            //TODO: Handle the exception appropriately
        }
    }
    
    // the onMessage method gets executed if a message arrives on the given channel...
	@Override
	public void onMessage(Message request) {
		//MessageProcessor.get(this.session, request).run();
		System.out.println("[ChannelListener.onMessage] -> "+request);
		new Thread(
				new MessageProcessor(this.session, request)).start();
	}	
}
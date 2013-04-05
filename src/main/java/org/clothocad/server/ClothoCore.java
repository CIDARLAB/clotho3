package org.clothocad.server;

import javax.jms.Connection;
import javax.jms.JMSException;
import javax.jms.Message;
import javax.jms.MessageConsumer;
import javax.jms.MessageListener;
import javax.jms.Session;

import org.apache.activemq.ActiveMQConnectionFactory;

public class ClothoCore 
		implements MessageListener {
	
	private Connection connection = null;
	
	public ClothoCore() 
			throws Exception {
		
		// first, we start the broker ...
		//new ClothoBroker().start();		
		
		// second, we connect to the broker ...
		this.connect();
		
		// third, we create a session for the connection
        Session session = connection.createSession(false, Session.CLIENT_ACKNOWLEDGE);

        // fourth, we create a consumer which reads incoming messages
        MessageConsumer consumer = session.createConsumer(
        		session.createQueue("CLOTHO"));
        
        // finally, we listen for incoming messages
        consumer.setMessageListener(this);  
        		// for every incoming message, 
        		// the onMessage method gets invoked      
    }
	
	/**
	protected ConnectionFactory createConnectionFactory() throws Exception {
        ConnectionFactory result = new ConnectionFactory();
        result.setBrokerURI("tcp://localhost:61613");
        return result;
    }
	**/
	
	private void connect() 
			throws Exception {		
		ActiveMQConnectionFactory cf = new ActiveMQConnectionFactory("admin", "password", "tcp://localhost:61613");
		this.connection = cf.createConnection();
		this.connection.start();
		//this.connection.open("localhost", 61613);

		/**
		this.connection = createConnectionFactory().createConnection();
        this.connection.setClientID("CLOTHO-CORE");
        this.connection.start();
        **/        
	}	

    // the onMessage method gets executed if a message arrives on the given channel...
	@Override
	public void onMessage(Message message) {
		//MessageProcessor.get(this.session, request).run();
		System.out.println("[ChannelListener.onMessage] -> "+message);
		
		/**
		StompJmsTextMessage txt = (StompJmsTextMessage)message;
		try {
			Enumeration<String> propNames = txt.getAllPropertyNames();
			while(propNames.hasMoreElements()) {
				System.out.println(propNames.nextElement());
			}
		} catch (JMSException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		**/
		
		//new Thread(
		//		new MessageProcessor(this.session, request)).start();
	}
	
	
	public void shutdown() 
			throws JMSException {
		if(null != connection) {
	        connection.close();
		}
	}
}

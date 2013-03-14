package org.clothocad.server;

import org.apache.activemq.broker.BrokerFactory;
import org.apache.activemq.broker.BrokerService;
import org.clothocad.core.layers.communication.ClothoConstants;
import org.clothocad.core.layers.communication.activemq.ChannelListener;

public class ClothoServer {
	private BrokerService broker = null;
	private ChannelListener listener = null;
	
	public ClothoServer() {
		this.broker = new BrokerService(); 
		/** todo: using a configuration file 
		this.broker = BrokerFactory.createBroker(
				"xbean:src/release/conf/activemq.xml");
		**/
		
		// here is the broker
        this.broker.setPersistent(false); // non-persistent
        this.broker.setUseJmx(true);        
    }
	
	public void start() {
    	try {
            this.broker.addConnector(ClothoConstants.SERVER_URL); // this is the server's URL                    
            this.broker.start(); // start the Broker

            
            System.out.println("Clotho's Broker is now running at "+ClothoConstants.SERVER_URL+"...");
            
           	this.listener = 
           			new ChannelListener(ClothoConstants.SERVER_URL);

    	} catch (Exception e) {
            e.printStackTrace();
        }    	
	}
	
	public void shutdown() {
		if(null != this.listener) {
			this.listener = (ChannelListener)null;
		}
		
		if (null != this.broker) {
			try {
				this.broker.stop();
			} catch (Exception e) {
				e.printStackTrace();
			}
			this.broker = (BrokerService)null;
		}
	}
}

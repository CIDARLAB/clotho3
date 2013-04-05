package org.clothocad.server;

import org.apache.activemq.broker.BrokerService;
import org.apache.activemq.broker.TransportConnector;

public class ClothoBroker {

	private BrokerService broker = null;
	
	public ClothoBroker() {
		this.broker = new BrokerService();
        this.broker.setPersistent(false);
	}
	
	public void start() 
			throws Exception {
        broker.addConnector("stomp://localhost:61613");
        broker.addConnector("tcp://0.0.0.0:61613");
        broker.start();
        broker.waitUntilStarted();
	}
}

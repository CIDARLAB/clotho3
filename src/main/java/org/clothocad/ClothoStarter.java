package org.clothocad;

import org.clothocad.broker.ClothoBroker;
import org.clothocad.core.ClothoCore;
import org.clothocad.webserver.jetty.ClothoWebserver;

public class ClothoStarter {
	public static void main(String[] args) 
			throws Exception {

		// start the message broker
		new ClothoBroker();
		
		// wait a bit until the broker is running
		Thread.sleep(4000);
		
		// now, start the core (i.e. the server)
		new ClothoCore().start();
		
		// wait a bit until the broker is running
		Thread.sleep(4000);
		
		// start the Jetty webserver
		new ClothoWebserver();
		
		Object lock = new Object();
        synchronized (lock) {
            lock.wait();
        }

	}
}

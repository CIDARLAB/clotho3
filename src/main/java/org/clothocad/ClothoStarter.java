package org.clothocad;

import org.clothocad.broker.ClothoBroker;
import org.clothocad.core.ClothoCore;

public class ClothoStarter {
	public static void main(String[] args) 
			throws Exception {

		new ClothoBroker();
		
		// wait a bit until the broker is running
		Thread.sleep(2000);
		
		new ClothoCore().start();	
		
		Object lock = new Object();
        synchronized (lock) {
            lock.wait();
        }

	}
}

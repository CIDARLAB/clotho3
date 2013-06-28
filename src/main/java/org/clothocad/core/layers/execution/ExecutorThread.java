package org.clothocad.core.layers.execution;

import java.util.Map;
import org.clothocad.core.layers.communication.Communicator;

public class ExecutorThread 
	implements Runnable {

	private String socket_id;
	private String channel;
	private Map<String, Object> json;
	
	public ExecutorThread(String socket_id, String channel, Map<String, Object> json) {
		this.socket_id = socket_id;
		this.channel = channel;
		this.json = json;
	}
	
	public void run() {
		// here we need to invoke the Executor...
		System.out.println("[ExecutorThread.run] -> "+socket_id);
		
		try {
			Thread.sleep(5000);
		} catch(Exception e) {
			// do nothing
		}
		
		// afterwards, respond to the client
		Communicator.get().sendClientMessage(socket_id, channel, "CORE-RESPONSE");
	}
}

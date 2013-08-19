package org.clothocad.core.testers.publisher;

import org.clothocad.core.layers.communication.apollo.ClothoPublisher;

import java.util.HashMap;
import java.util.Map;

public class SendMessagesToClient {
	
	public SendMessagesToClient() {
		for(int i=1; i<=10; i++) {
			Map<String, Object> json = new HashMap<>();
			try {
				json.put("update", "UPDATE");
	
				// send messages to the CLOTHO.UPDATES topic
				new ClothoPublisher().publish(json);
			} catch (Exception e) {
				e.printStackTrace();
			}
			
			try {
				Thread.sleep(5000);
			} catch(Exception e) {}
		}
	}
	
	public static void main(String[] args) {
		new SendMessagesToClient();
	}
	
}

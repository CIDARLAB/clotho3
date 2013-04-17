package org.clothocad.core.testers.publisher;

import org.clothocad.core.layers.communication.activemq.ClothoPublisher;
import org.json.JSONException;
import org.json.JSONObject;

public class SendMessagesToClient {
	
	public SendMessagesToClient() {
		for(int i=1; i<=10; i++) {
			JSONObject json = new JSONObject();
			try {
				json.put("update", "UPDATE");
	
				// send messages to the CLOTHO.UPDATES topic
				new ClothoPublisher().publish(json);
			} catch (JSONException e) {
				e.printStackTrace();
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

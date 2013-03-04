package org.clothocad.core.layers.communication;

import org.json.JSONObject;

public class Publisher {
	
	private static Publisher publisher = (Publisher)null;
	
	public static Publisher get() {
		if((Publisher)null == publisher) {
			publisher = new Publisher();
		}
		
		return publisher;
	}

	public void publish(JSONObject json) {
		try {
			Router.get().sendMessage("", Channel.UPDATES.toString(), json);
		} catch(Exception e) {
			e.printStackTrace();
		}
	}
	
}

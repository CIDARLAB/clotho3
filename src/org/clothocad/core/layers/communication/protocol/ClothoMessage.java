package org.clothocad.core.layers.communication.protocol;

import org.json.JSONObject;

public class ClothoMessage {
	private JSONObject requestData;
	
	// the JSON object must consist of...
	// 1. the channel name (i.e. CHANNEL key)
	// 2. the action name (i.e. ACTION key)
	// 3. a datum object (i.e. DATUM key)
	// 4. the client's IP address (i.e. IP key)
	// 5. the user's id (i.e. the USER key)
	
	public ClothoMessage(JSONObject requestData) {
		this.requestData = requestData;
	}
	
	public JSONObject getRequestData() {
		return this.requestData;
	}
}

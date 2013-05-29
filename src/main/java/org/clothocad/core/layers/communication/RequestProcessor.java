package org.clothocad.core.layers.communication;

import org.clothocad.core.layers.communication.connection.ClientConnection;
import org.clothocad.core.layers.communication.connection.ws.ClothoWebSocket;
import org.clothocad.core.layers.communication.protocol.ActionType;
import org.json.JSONException;
import org.json.JSONObject;

import ch.qos.logback.core.joran.action.Action;

public class RequestProcessor 
	implements Runnable {

	private ClientConnection connection;
	private String action;
	private String correlationId;
	private JSONObject data;
	
	public RequestProcessor(ClientConnection connection, String action, String correlationId, JSONObject data) {
		this.connection = connection;
		this.action = action;
		this.correlationId = correlationId;
		this.data = data;
	}
	
	@Override
	public void run() {
	
		try {
			if(ActionType.GET.toString().equalsIgnoreCase(action)) {
				// do the request processing
				new ServerSideAPI(connection).get(this.data.getString("uuid"));
			} else if(ActionType.SET.toString().equalsIgnoreCase(action)) {
				// do the request processing
				//new ServerSideAPI(connection).set(data);
			} else if(ActionType.SAY.toString().equalsIgnoreCase(action)) {
				
			} else if(ActionType.SUBMIT.toString().equalsIgnoreCase(action)) {
				
			} else if(ActionType.RUN.toString().equalsIgnoreCase(action)) {
				
			} else if(ActionType.CREATE.toString().equalsIgnoreCase(action)) {
				
				String uuid = new ServerSideAPI().create(this.data);
				
				// create the JSONObject to be responded
				JSONObject response = new JSONObject();

				// header-information
				response.put(ClothoConstants.CHANNEL, Channel.RESPONSE);
				response.put(ClothoConstants.CORRELATION_ID, correlationId);
				
				// data
				JSONObject responseData = new JSONObject();
				responseData.put("uuid", uuid);
				response.put(ClothoConstants.DATA, responseData);
				
				// return the response to the requesting client
				Router.get().sendMessage(connection, Channel.RESPONSE.toString(), response);
			}
		} catch (JSONException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}

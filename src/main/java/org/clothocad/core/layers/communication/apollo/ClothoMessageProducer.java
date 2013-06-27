package org.clothocad.core.layers.communication.apollo;

import org.clothocad.core.layers.communication.Callback;
import org.clothocad.core.layers.communication.connection.apollo.ApolloConnection;
import org.json.JSONObject;

public class ClothoMessageProducer 
		implements Callback {
	
	private ApolloConnection connection;
	public ClothoMessageProducer(ApolloConnection connection) {
		this.connection = connection;
	}

	@Override
	public void onSuccess(JSONObject json) {
            CallbackHandler cbh = 
                    CallbackHandlerTable.get(connection.getId());
            if(null != cbh) {
                System.err.println("[ClothoMessageProducer.onSuccess] -> "+json+" -> "+connection.getId());
                cbh.respond(json);			
            }
	}

	@Override
	public void onFailure(Throwable err) {
            CallbackHandler cbh = CallbackHandlerTable.get(
                            connection.getId());
            if(null != cbh) {
                JSONObject jsonResponse = new JSONObject();
                try {
                    jsonResponse.put("CLOTHO-ERROR", err.getLocalizedMessage());
                } catch(Exception e) {}
                System.err.println("[ClothoMessageProducer.onFailure] -> "+jsonResponse+" -> "+connection.getId());
                cbh.respond(jsonResponse);			
            }
	}

}

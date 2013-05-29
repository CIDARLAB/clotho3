package org.clothocad.core.layers.communication.activemq;

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
		CallbackHandler cbh = CallbackHandlerTable.get(
				connection.getCorrelationId());
		if(null != cbh) {
			cbh.respond(json);			
		}
	}

	@Override
	public void onFailure(Throwable err) {
		// TODO Auto-generated method stub
		
	}

}

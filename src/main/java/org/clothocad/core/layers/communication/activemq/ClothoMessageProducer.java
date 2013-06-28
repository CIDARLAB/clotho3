package org.clothocad.core.layers.communication.activemq;

import java.util.Map;
import org.clothocad.core.layers.communication.Callback;
import org.clothocad.core.layers.communication.connection.apollo.ApolloConnection;

public class ClothoMessageProducer 
		implements Callback {
	
	private ApolloConnection connection;
	public ClothoMessageProducer(ApolloConnection connection) {
		this.connection = connection;
	}

	@Override
	public void onSuccess(Object json) {
		CallbackHandler cbh = CallbackHandlerTable.get(
				connection.getCorrelationId());
		if(null != cbh) {
			System.err.println("[ClothoMessageProducer.onSuccess] -> "+json);
			cbh.respond(json);			
		}
	}

	@Override
	public void onFailure(Throwable err) {
		// TODO Auto-generated method stub
		
	}

}

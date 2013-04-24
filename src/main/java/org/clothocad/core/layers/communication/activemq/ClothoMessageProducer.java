package org.clothocad.core.layers.communication.activemq;

import org.clothocad.core.layers.communication.Callback;
import org.json.JSONObject;

public class ClothoMessageProducer 
		implements Callback {
	
	private String sCorrelationID;
	
	public ClothoMessageProducer(String sCorrelationID) {
		this.sCorrelationID = sCorrelationID;
	}

	@Override
	public void onSuccess(JSONObject json) {
		CallbackHandler cbh = CallbackHandlerTable.get(sCorrelationID);
		if(null != cbh) {
			cbh.respond(json);			
		}
	}

	@Override
	public void onFailure(Throwable err) {
		// TODO Auto-generated method stub
		
	}

}

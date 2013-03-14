package org.clothocad.core.layers.communication.activemq;

import org.clothocad.core.layers.communication.Callback;
import org.json.JSONObject;

public class ClothoMessageProducer 
		implements Callback {
	
	private String sCallbackHandlerID;
	
	public ClothoMessageProducer(String sCallbackHandlerID) {
		this.sCallbackHandlerID = sCallbackHandlerID;
	}

	@Override
	public void onSuccess(JSONObject json) {
		CallbackHandler cbh = CallbackHandlerTable.get(sCallbackHandlerID);
		if(null != cbh) {
			cbh.respond(json);
			
			// then, remove the Callback handler from the Callback-handler table...
			CallbackHandlerTable.remove(cbh);
		}
	}

	@Override
	public void onFailure(Throwable err) {
		// TODO Auto-generated method stub
		
	}
}

package org.clothocad.core.layers.communication.activemq;

import java.util.Hashtable;

public class CallbackHandlerTable {

	// key   ... the message's correlation id
	// value ... the callback-handler object 
	private static Hashtable<String, CallbackHandler> htCallbackHandlers;
	
	public static void put(String sCorrelationID, CallbackHandler cbh) {
		if(null == htCallbackHandlers) {
			htCallbackHandlers = new Hashtable<String, CallbackHandler>();
		}
		
		if(!htCallbackHandlers.containsKey(sCorrelationID)) {
			htCallbackHandlers.put(sCorrelationID, cbh);
		}
	}
	
	public static CallbackHandler get(String sCorrelationID) {
		if(null != htCallbackHandlers) {
			CallbackHandler cbh = htCallbackHandlers.get(sCorrelationID);
			htCallbackHandlers.remove(sCorrelationID);
			return cbh;
		}
		return (CallbackHandler)null;
	}
}

package org.clothocad.core.layers.communication.activemq;

import java.util.Hashtable;

public class CallbackHandlerTable {

	private static Hashtable<String, CallbackHandler> htSessions;
	
	public static String put(CallbackHandler cbh) {
		if(null == htSessions) {
			htSessions = new Hashtable<String, CallbackHandler>();
		}
		
		String sKey = String.valueOf(cbh.hashCode());
		if(!htSessions.containsKey(sKey)) {
			htSessions.put(sKey, cbh);
		}
		return sKey;
	}
	
	public static CallbackHandler get(String sKey) {
		if(null != htSessions) {
			return htSessions.get(sKey);
		}
		return (CallbackHandler)null;
	}
	
	public static void remove(CallbackHandler cbh) {
		if (null != htSessions) {
			String sKey = String.valueOf(cbh.hashCode());
			if(htSessions.containsKey(sKey)) {
				htSessions.remove(sKey);
			}
		}
	}
}

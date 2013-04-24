package org.clothocad.core.layers.communication.protocol;

import java.net.InetAddress;

import org.clothocad.core.datums.Datum;
import org.json.JSONObject;

public class Marshaller {

	public static JSONObject marshal(ActionType action, JSONObject json) {
		JSONObject ret = new JSONObject();
		try {
			ret.put("ip", InetAddress.getLocalHost());
			ret.put("action", action.toString());
			ret.put("datum", json);
		} catch (Exception e) {
			e.printStackTrace();
			return (JSONObject)null;
		}
		return json;
	}

	public static JSONObject toJSON(ActionType action, Datum datum) {
		JSONObject json = new JSONObject();
		try {
			json.put("ip", InetAddress.getLocalHost());
			json.put("action", action.toString());
			
			// TODO: serialize the datum to JSON
			//json.put("datum", new JSONSerializer().t datum.toJSON());
		} catch (Exception e) {
			e.printStackTrace();
			return (JSONObject)null;
		}
		return json;
	}
}

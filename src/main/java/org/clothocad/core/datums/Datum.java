package org.clothocad.core.datums;

import java.io.Serializable;
import java.util.UUID;

import org.clothocad.core.datums.util.ClothoDate;
import org.json.JSONException;
import org.json.JSONObject;

public interface Datum 
	extends Serializable {
	
	String getId();
	
}

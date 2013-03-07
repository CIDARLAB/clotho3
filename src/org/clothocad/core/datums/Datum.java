package org.clothocad.core.datums;

import java.io.Serializable;
import java.util.UUID;

import org.clothocad.core.datums.util.ClothoDate;
import org.json.JSONException;
import org.json.JSONObject;

public abstract class Datum 
	implements Serializable {
	
	private static final long serialVersionUID = -4383410157216396194L;
	
	private UUID uuid;
	
	private ClothoDate dateCreated;
	private ClothoDate lastModified;
	private ClothoDate lastAccessed;
	
	
	protected JSONObject json;
	public Datum(String json) 
			throws JSONException {
		this.uuid = UUID.randomUUID();
		this.dateCreated = new ClothoDate();
		this.lastModified = new ClothoDate();
		this.lastAccessed = new ClothoDate();
		this.json = new JSONObject(json);
	}
	
	public Datum() {
		// create a random UUID
		this.uuid = UUID.randomUUID();
		this.dateCreated = new ClothoDate();
		this.lastModified = new ClothoDate();
		this.lastAccessed = new ClothoDate();
		this.json = new JSONObject();
	}	
	
	/**
	public UUID getUUID() {
		return this.uuid;
	}
	**/
	
	public String getId() {
		return this.uuid.toString();
	}
	
	public ClothoDate getDateCreated() {
		return this.dateCreated;
	}
	
	public void setLastModified(ClothoDate lastModified) {
		this.lastModified = lastModified;
	}
	
	public ClothoDate getLastModified() {
		return this.lastModified;
	}
	
	
	public void setLastAccessed(ClothoDate lastAccessed) {
		this.lastAccessed = lastAccessed;
	}
	
	public ClothoDate getLastAccessed() {
		return this.lastAccessed;
	}
	
    // every Datum must implement the toJSON method
    public JSONObject toJSON()  {   
    	if (null != json) {
    		return this.json;
    	}
    	// here we can utilize Stephanie's introspection code to 
    	// map datums into JSON/BSON format
    	return new JSONObject();
    }
}

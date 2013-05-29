/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.datums;

import org.clothocad.model.Person;
import org.json.JSONObject;

/**
 *
 * @author jcanderson
 */
public abstract class Sharable 
	extends ObjBase 
	implements Datum {
    
	private Person author;
	
	public Sharable(String name, Person author) {
		super(name);
		this.author = author;
	}
	
	public Sharable(String json) 
			throws Exception {
		new JSONObject(json);
		// TODO: parse the JSON and set the appropriate data fields (i.e. author)
	}
	
	public Person getAuthor() {
		return this.author;
	}

	public abstract JSONObject toJSON();
    public abstract Person extractAuthor();
    //NEED BACK THE SHARING RULES FROM OLD DATUM
    
    public abstract SharableType type();

    /**
     * Each Sharable class is responsible for handling its own permissions
     * 
     * @param newvalue
     * @param requestor
     * @return 
     */
    public abstract boolean set(JSONObject newvalue, Person requestor, Doo doo);

    public static enum SharableType {SCHEMA, VIEW, FUNCTION, INSTANCE, BADGE, QUIZ, TRAIL, CONTENT};
}

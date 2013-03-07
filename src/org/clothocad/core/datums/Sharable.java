/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.datums;

import org.clothocad.core.aspects.Collector;
import org.clothocad.core.datums.objbases.Person;

/**
 *
 * @author jcanderson
 */
public abstract class Sharable 
		extends Datum {
	
	private Person author;
	private SharableType type;
	
	public Sharable() {}
	
	public Sharable(Person author, SharableType type) {
		super();
		this.author = author;
		this.type = type;
	}
	
	
    public Person extractAuthor() {
    	return this.author;
    	
    	// 2nd option:
    	//     return (Person) Collector.get().getDatum(authorId);
    }
    
    public SharableType getType() {
    	return this.type;
    }

    
    // ???
    //public abstract boolean set(JSONObject newvalue, Person requestor, Doo doo);
    
    /***
    JSONObject toJSON();
    //NEED BACK THE SHARING RULES FROM OLD DATUM
    
    public SharableType type();

    /**
     * Each Sharable class is responsible for handling its own permissions
     * 
     * @param newvalue
     * @param requestor
     * @return 
     */
	/***
    public boolean set(JSONObject newvalue, Person requestor, Doo doo);
    ***/

    public static enum SharableType {SCHEMA, VIEW, FUNCTION, INSTANCE, BADGE, QUIZ, TRAIL, CONTENT};
}

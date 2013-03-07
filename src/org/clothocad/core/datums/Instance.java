/*
 * 
Copyright (c) 2010 The Regents of the University of California.
All rights reserved.
Permission is hereby granted, without written agreement and without
license or royalty fees, to use, copy, modify, and distribute this
software and its documentation for any purpose, provided that the above
copyright notice and the following two paragraphs appear in all copies
of this software.

IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.

THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
ENHANCEMENTS, OR MODIFICATIONS..
 */

package org.clothocad.core.datums;

import org.clothocad.core.aspects.Collector;
import org.clothocad.core.datums.objbases.Person;
import org.clothocad.core.datums.util.Permissions;
import org.clothocad.core.util.Logger;
import org.json.JSONException;
import org.json.JSONObject;

/**
 * @author John Christopher Anderson
 */


public class Instance 
	extends Sharable {

	private static final long serialVersionUID = 7626966374337532946L;

	private Schema schema;
    private Permissions perms = new Permissions();

	private Instance(Person author, Schema schema, String json) 
			throws JSONException {
		super(author, SharableType.INSTANCE);
		this.json = new JSONObject(json);
		
		this.schema = schema;
	}


    public static Instance create(Person author, Schema schema, String json) 
    		throws JSONException {
    	Instance inst = new Instance(author, schema, json);
    	Collector.get().add(inst);
    	return inst;
    }


    public Schema getSchema() {
    	return this.schema;
    }
    
    
    public String getString(String key) {
    	if(null != this.json && null != key) {
    		try {
				return json.getString(key);
			} catch (JSONException e) {
	            Logger.log(Logger.Level.WARN, key+" not found!", e);
			}
    	}
    	return null;
    }
    
    public static Instance deserialize(String jsonstr) {
        try {
            return new Instance(null, null, jsonstr);
        } catch (Exception ex) {
            return null;
        }
    }
    
    /**
    @Override
    public String getId() {
        try {
            return super.getString("id");
        } catch(Exception err) {
            return null;
        }
    }
    **/
    
    /***
    @Override
    public boolean set(JSONObject newvalue, Person requestor, Doo doo) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    ***/    
}
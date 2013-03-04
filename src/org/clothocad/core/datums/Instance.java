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

import org.clothocad.core.datums.util.ClothoDate;
import org.clothocad.core.datums.util.Permissions;
import org.clothocad.core.datums.objbases.Person;
import java.util.UUID;
import org.clothocad.core.aspects.Collector;
import org.codehaus.jettison.json.JSONException;
import org.json.JSONObject;

/**
 * @author John Christopher Anderson
 */


public class Instance 
	extends JSONObject 
	implements Sharable {

	private Schema schema;
	
    public Instance() {}
    
    protected Instance(String json) 
    		throws JSONException, org.json.JSONException {
        super(json);
    }
    
    public Instance(String json, Schema schema) 
    		throws Exception {
        super(json);
        this.schema = schema;
    }

    public static Instance create(Person author, Schema schema, String json) {
        try {
        	
            Instance out = new Instance(json, schema);
            out.put("schemaId", schema.getId());
            if(author!=null) {
                out.put("authorId", author.getId());
            }
            out.put("id", UUID.randomUUID().toString());
            ClothoDate date = new ClothoDate();
            out.put("dateCreated", date.toJSON());
            out.put("dateLastModified", date.toJSON());
            out.put("dateLastAccessed", date.toJSON());
            Permissions perms = new Permissions();
            out.put("permissions", perms.toJSON());
            return out;
        } catch(Exception err) {
            err.printStackTrace();
            return null;
        }
    }


    /**
    public Schema getSchema() {
    	return this.schema;
    }
    **/
    
    @Override
    public Person extractAuthor() {
        try {
            String authorId = super.getString("authorId");
            Person out = (Person) Collector.get().getDatum(authorId);
            return out;
        } catch(Exception err) {
            err.printStackTrace();
            return null;
        }
    }
    
    @Override
    public JSONObject toJSON() {
        JSONObject out = new JSONObject();
        return this;
    }

    public static Instance deserialize(String jsonstr) {
        try {
            return new Instance(jsonstr);
        } catch (Exception ex) {
            return null;
        }
    }
    
    @Override
    public String getId() {
        try {
            return super.getString("id");
        } catch(Exception err) {
            return null;
        }
    }

    @Override
    public boolean set(JSONObject newvalue, Person requestor, Doo doo) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public SharableType type() {
        return SharableType.INSTANCE;
    }
}
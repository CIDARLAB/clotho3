/*
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
ENHANCEMENTS, OR MODIFICATIONS.
 */

package org.clothocad.core.datums;

import org.clothocad.core.datums.util.ServerScript;
import org.clothocad.core.datums.util.ClothoDate;
import org.clothocad.core.datums.util.Permissions;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.datums.objbases.Person;
import flexjson.JSONDeserializer;
import flexjson.JSONSerializer;
import java.util.List;
import java.util.Map;
import java.util.UUID;
import org.clothocad.core.aspects.Collector;
import org.json.JSONObject;



/**
 * @author John Christopher Anderson
 */

/**
 * A Schema is the schema that must be obeyed by an Instance
 * @author jcanderson
 */
public class Schema 
		extends Sharable {

	public Schema() {}
	
    private Schema(Person author,
            String name, 
            String description,
            List<ClothoField> fields,
            ServerScript indexer,
            ServerScript queryer,
            ServerScript validator) {
    	
    	super(author, SharableType.SCHEMA);
    	
    	this.name = name;
    	this.fields = fields;    	
        this.name = name;
        this.description = description;
        this.fields = fields;
        this.indexingScript = indexer;
        this.validationScript = validator;
        this.queryingScript = queryer;

    }
    
    /**
     * Instantiating a new Schema from it's various components
     * @param name
     * @param fields
     * @param validators
     * @param indexer  //This script provides the logic needed to customize the collator somehow, or does it's own custom indexing
     * @param example
     * @param prototype
     * @return 
     */
    public static Schema create(Person author,
                                String name, 
                                String description,
                                List<ClothoField> fields,
                                ServerScript indexer,
                                ServerScript queryer,
                                ServerScript validator) {
        
        //CHECK THE DATA FOR WELL-FORMEDNESS
        
        //Return the schema
        Schema schema = new Schema(author, name, description, fields, indexer, queryer, validator);
        
        Collector.get().add(schema);
        return schema;
    }

    public boolean validate(String jsondata) {
        return false;
    }

    @Override
    public Person extractAuthor() {
        Person out = (Person) Collector.get().getDatum(authorId);
        return out;
    }
    
    /***
    @Override
    public JSONObject toJSON() {
        try {
            JSONSerializer serializer = new JSONSerializer().exclude("*.class");
            serializer.prettyPrint(true);
            String serial = serializer.deepSerialize( this );
            return new JSONObject(serial);
        } catch (Exception ex) {
            return null;
        }
    }
    ***/
    
    public static Schema deserialize(String json) {
        return new JSONDeserializer<Schema>().deserialize(json, Schema.class);
    }

    public String getDescription() {
        return description;
    }

    public List<ClothoField> getFields() {
        return fields;
    }

    public ServerScript getIndexingScript() {
        return indexingScript;
    }

    public String getName() {
        return name;
    }

    public ServerScript getQueryingScript() {
        return queryingScript;
    }

    public ServerScript getValidationScript() {
        return validationScript;
    }

    public String getAuthorId() {
        return authorId;
    }

    public int getInstanceCount() {
        return instanceCount;
    }

    public String getLargeIconURL() {
        return largeIconURL;
    }

    public Permissions getPermissions() {
        return permissions;
    }

    public String getSmallIconURL() {
        return smallIconURL;
    }

    public String getViewId() {
        return viewId;
    }

    /***
    @Override
    public boolean set(JSONObject newvalue, Person requestor, Doo doo) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    ***/
    
    private ServerScript indexingScript;
    private ServerScript queryingScript;
    private ServerScript validationScript;
    
    private List<ClothoField> fields;
    
    //Permissions
    private Permissions permissions = new Permissions();

    
    //Metadata
    private String name;
    private String description;
    
    private String viewId;
    private String authorId;
    private String smallIconURL;
    private String largeIconURL;
    
    private int instanceCount = 0;
}


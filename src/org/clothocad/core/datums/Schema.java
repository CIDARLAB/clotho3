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
 * A Schema is the schema that must be obeyed by an ObjBase
 * @author jcanderson
 */
public class Schema implements Sharable {

    private Schema() {
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
        Schema out = new Schema();

        if(author!=null) {   //PROBABLY SHOULD GET RID OF THIS LATER
            out.authorId = author.getId();
        }
        out.name = name;
        out.id = UUID.randomUUID().toString();
        out.description = description;
        out.fields = fields;
        out.indexingScript = indexer;
        out.validationScript = validator;
        out.queryingScript = queryer;

        Collector.get().add(out);
        return out;
    }

    public boolean validate(String jsondata) {
        return false;
    }

    @Override
    public Person extractAuthor() {
        Person out = (Person) Collector.get().getDatum(authorId);
        return out;
    }
    
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
    
    public static Schema deserialize(String json) {
        Schema out = new JSONDeserializer<Schema>().deserialize(json, Schema.class);
        return out;
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

    public ClothoDate getDateCreated() {
        return dateCreated;
    }

    public ClothoDate getDateLastAccessed() {
        return dateLastAccessed;
    }

    public ClothoDate getDateLastModified() {
        return dateLastModified;
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
    
    @Override
    public String getId() {
        return id;
    }

    @Override
    public SharableType type() {
        return SharableType.SCHEMA;
    }

    @Override
    public boolean set(JSONObject newvalue, Person requestor, Doo doo) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    private ServerScript indexingScript;
    private ServerScript queryingScript;
    private ServerScript validationScript;
    
    private List<ClothoField> fields;
    
    //Permissions
    private Permissions permissions = new Permissions();

    
    //Metadata
    private String id;
    private String name;
    private String description;
    
    private String viewId;
    private String authorId;
    private String smallIconURL;
    private String largeIconURL;
    
    private int instanceCount = 0;
    
    private ClothoDate dateCreated = new ClothoDate();
    private ClothoDate dateLastModified = new ClothoDate();
    private ClothoDate dateLastAccessed = new ClothoDate();
}


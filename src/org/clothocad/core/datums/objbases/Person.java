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

package org.clothocad.core.datums.objbases;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.clothocad.core.aspects.Collector;
import org.clothocad.core.datums.Datum;
import org.clothocad.core.datums.Instance;
import org.clothocad.core.datums.Schema;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.datums.util.FieldType;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.datums.util.ServerScript;
import org.json.JSONObject;

/**
 * @author John Christopher Anderson
 */
public class Person 
		extends Datum {

	private static final long serialVersionUID = -5721110846830687539L;

	private Person(String json) 
    		throws Exception {
        super(json);
    }
	
    /**
    public String getName() {
        try {
            return this.getString("displayname");
        } catch (org.json.JSONException ex) {
            return "exception";
        }
    }
    **/
    
	/**
    public static Person create(Instance obj) {
        try {
            JSONObject json = obj.toJSON();
            String jsonstr = json.toString();
            
            Person out = new Person(jsonstr);
            return out;
        } catch (Exception ex) {
            return null;
        }
    }
    **/
	
    public static Person getById(String uuid) {
        try {
            Instance obj = (Instance) Collector.get().getDatum(uuid);
            JSONObject json = obj.toJSON();
            String jsonstr = json.toString();
            
            Person out = new Person(jsonstr);
            return out;
        } catch (Exception ex) {
            return null;
        }
    }	
    
    public static Schema getSchema() {
    	if(personSchema == null) {
            String json = ""
                    + "{\"authorId\":\"static-admin-instance-is-uuid\",\"dateCreated\":{\"absolute\":1345997642813,\"day\":26,\"hour\":9,\"millis\":813,\"minute\":14,\"month\":\"August\",\"second\":2,\"year\":2012},\"dateLastAccessed\":{\"absolute\":1345997642814,\"day\":26,\"hour\":9,\"millis\":814,\"minute\":14,\"month\":\"August\",\"second\":2,\"year\":2012},\"dateLastModified\":{\"absolute\":1345997642814,\"day\":26,\"hour\":9,\"millis\":814,\"minute\":14,\"month\":\"August\",\"second\":2,\"year\":2012},\"description\":\"The Sharable representation of a singular human being.\",\"fields\":{\"middlename\":{\"example\":\"Marcus\",\"permissions\":1,\"type\":\"STRING\"},\"nickname\":{\"example\":\"Bob\",\"permissions\":1,\"type\":\"STRING\"},\"email\":{\"example\":\"bob@bob.com\",\"permissions\":1,\"type\":\"STRING\"},\"givenname\":{\"example\":\"William\",\"permissions\":1,\"type\":\"STRING\"},\"surname\":{\"example\":\"Fiddleston\",\"permissions\":1,\"type\":\"STRING\"},\"displayname\":{\"example\":\"bob.fiddleston\",\"permissions\":1,\"type\":\"STRING\"}},\"id\":\"static-person-schema-is-uuid\",\"indexingScript\":{\"language\":\"JavaScript\",\"script\":\"return true;\"},\"instanceCount\":0,\"largeIconURL\":null,\"name\":\"Person\",\"permissions\":{\"entries\":[],\"restOfWorld\":0},\"queryingScript\":{\"language\":\"JavaScript\",\"script\":\"return true;\"},\"smallIconURL\":null,\"validationScript\":{\"language\":\"JavaScript\",\"script\":\"return true;\"},\"viewId\":null}"
                    + "";
            personSchema = Schema.deserialize(json);
        }
        return personSchema;
    }    

    public static Person getAdmin() {
        try {
            if(adminPerson == null) {
                adminPerson = new Person(""
                        + "{\"email\":\"yourhandle@yourdomain.com\",\"displayname\":\"admin\",\"givenname\":\"Admin\",\"middlename\":\"I\",\"surname\":\"Strator\",\"nickname\":\"admin\",\"schemaId\":\"static-person-schema-is-uuid\",\"id\":\"static-admin-instance-is-uuid\",\"dateCreated\":{\"absolute\":1345997642929,\"day\":26,\"hour\":9,\"millis\":929,\"minute\":14,\"month\":\"August\",\"second\":2,\"year\":2012},\"dateLastModified\":{\"absolute\":1345997642929,\"day\":26,\"hour\":9,\"millis\":929,\"minute\":14,\"month\":\"August\",\"second\":2,\"year\":2012},\"dateLastAccessed\":{\"absolute\":1345997642929,\"day\":26,\"hour\":9,\"millis\":929,\"minute\":14,\"month\":\"August\",\"second\":2,\"year\":2012},\"permissions\":{\"entries\":[],\"restOfWorld\":0}}"
                        + "");
            }            
            return adminPerson;
        } catch (Exception ex) {
            Logger.getLogger(Person.class.getName()).log(Level.SEVERE, null, ex);
            return null;
        }
    }
    
    public static Schema createStaticPersonClass() {
        List<ClothoField> fields = new ArrayList<ClothoField>();

        fields.add(new ClothoField("email", FieldType.STRING, "bob@bob.com", 1) );
        fields.add(new ClothoField("displayname", FieldType.STRING, "bob.fiddleston", 1) );
        fields.add(new ClothoField("givenname", FieldType.STRING, "William", 1) );
        fields.add(new ClothoField("middlename", FieldType.STRING, "Marcus", 1) );
        fields.add(new ClothoField("surname", FieldType.STRING, "Fiddleston", 1) );
        fields.add(new ClothoField("nickname", FieldType.STRING, "Bob", 1) );

        ServerScript script = new ServerScript("return true;", Language.JavaScript);
        Schema schema = Schema.create(null,
                "Person",
                "The Sharable representation of a singular human being.",
                fields,
                script,
                script,
                script);

        System.out.println(schema.toJSON().toString());

        return schema;
    }

    public static void createAdminPersonInstance() {
        try {
            JSONObject fields = new JSONObject();

            fields.put("email",  "yourhandle@yourdomain.com" );
            fields.put("displayname",  "admin" );
            fields.put("givenname",  "Admin" );
            fields.put("middlename",  "I" );
            fields.put("surname",  "Strator" );
            fields.put("nickname",  "admin" );

            String str = fields.toString();
            Instance obj = Instance.create(null, Person.createStaticPersonClass(), str);

            System.out.println(obj.toString());

        } catch (Exception ex) {
            Logger.getLogger(Person.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    /**
    public static void createRandomPersonInstance() {
        try {
            JSONObject fields = new JSONObject();

            fields.put("email",  "yourhandle@yourdomain.com" );
            fields.put("displayname",  "admin" );
            fields.put("givenname",  "Admin" );
            fields.put("middlename",  "I" );
            fields.put("surname",  "Strator" );
            fields.put("nickname",  "admin" );

            String str = fields.toString();
            Instance obj = Instance.create(null, Person.createStaticPersonClass(), str);

            System.out.println(obj.toString());

        } catch (Exception ex) {
            Logger.getLogger(Person.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
	**/
    
    private static Schema personSchema = null;
    private static Person adminPerson = null;
    
	@Override
	public JSONObject toJSON() {
		// TODO Auto-generated method stub
		return null;
	}
}
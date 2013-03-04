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
package org.clothocad.core.testers;

import flexjson.JSONSerializer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.clothocad.core.aspects.Collector;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.datums.util.FieldType;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.datums.Instance;
import org.clothocad.core.datums.objbases.Person;
import org.clothocad.core.datums.Schema;
import org.clothocad.core.datums.util.ServerScript;
import org.codehaus.jettison.json.JSONException;
import org.json.JSONObject;

/**
 * @author John Christopher Anderson
 */


public class Ta {
    public static void main(String[] args) throws Exception {
        //Create a schema and test
        System.out.println("\n\n+++++++++++ Testing out schema");
        Schema schema = createSchema();
        testSchema(schema);
        
        //Create and test an ObjBase
        System.out.println("\n\n+++++++++++ Testing out objbase");
        String serial = createFeature();
        Instance obj = Instance.create(Person.getAdmin(), schema, serial);
        System.out.println(obj.toJSON().toString());
        System.out.println("Id is: " + obj.getId());
        testObjBase(obj);
        
        //Try to load a Person object
        Person aperson = Person.getAdmin();
        System.out.println(aperson.toJSON().toString());
    }

    private static String createFeature() {
        HashMap afeat = new HashMap();
            afeat.put("name", "mRFP1");
            afeat.put("sequence", "atgcatcatcatcatcatcatta");
        
        JSONSerializer serializer = new JSONSerializer().exclude("*.class");
        serializer.prettyPrint(true);
        String serial = serializer.deepSerialize( afeat );
        System.out.println(serial);
        return serial;
    }
   
    
    private static Schema createSchema() throws JSONException {
        ServerScript indexer = new ServerScript("return true;", Language.JavaScript);
        ServerScript queryer = new ServerScript("return true;", Language.JavaScript);
        ServerScript validator = new ServerScript("return true;", Language.JavaScript);
        
        List<ClothoField> thefields = new ArrayList<ClothoField>();
        
        ClothoField name = new ClothoField("name", FieldType.NAME, "GFPmut3", 0);
        thefields.add(name);
        
        ClothoField sequence = new ClothoField("sequence", FieldType.STRING, "atgcatgagatcatgcagccaactatttattaa", 0);
        thefields.add(sequence);

        Schema schema = Schema.create(Person.getAdmin(), "Feature", "Act ontology standard representation of a genetic feature", thefields, indexer, queryer, validator);
        System.out.println(schema.toJSON().toString());
        return schema;
    }
    
    private static boolean testSchema(Schema schema) {
        try {
            //Convert to JSON
            JSONObject json = schema.toJSON();
            //Convert back to Schema
            Schema converted = Schema.deserialize(json.toString());

            System.out.println(converted.toJSON()); 
            
            
            if(!converted.getId().equals(schema.getId())) {
                return false;
            }
            
            
            System.out.println("Schema looks to be fine!!!");
            return true;
        } catch(Exception err) {
            err.printStackTrace();
            return false;  
        }
    }

    private static boolean testObjBase(Instance obj) {
        try {
            //Convert to JSON
            JSONObject json = obj.toJSON();
            //Convert back to Schema
            Instance converted = Instance.deserialize(json.toString());

            System.out.println(converted.toJSON()); 
            
            
            if(!converted.getId().equals(obj.getId())) {
                return false;
            }
            
            System.out.println("ObjBase looks to be fine!!!");
            return true;
        } catch(Exception err) {
            err.printStackTrace();
            return false;  
        }
    }
}

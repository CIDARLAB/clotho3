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
package org.clothocad.core.testers.schemas;

import static com.google.common.collect.Sets.newHashSet;
import com.mongodb.BasicDBObject;
import java.util.Set;
import org.bson.BSONObject;
import org.bson.types.ObjectId;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.datums.ObjBase;

import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.persistence.DBClassLoader;
import org.clothocad.core.persistence.mongodb.MongoDBConnection;
import org.clothocad.core.schema.Access;
import org.clothocad.core.schema.ClothoSchema;
import org.clothocad.core.schema.JavaSchema;
import org.clothocad.core.schema.Schema;
import static org.objectweb.asm.Opcodes.*;

/**
 * @author John Christopher Anderson
 */


public class SchemaATest {
    public static void main(String[] args) throws Exception {
        Persistor persistor = new Persistor(new MongoDBConnection());
        DBClassLoader cl = new DBClassLoader(persistor);
        
        //Create a schema and test
        System.out.println("\n\n+++++++++++ Testing out schema");
        Schema schema = createSchema();
        testSchema(schema);
        
        //Create and test an ObjBase
        System.out.println("\n\n+++++++++++ Testing out objbase");
        BSONObject serial = createFeature();
        persistor.save(serial.toMap());
        ObjBase feature = persistor.get(schema.getEnclosedClass(cl), new ObjectId((String) serial.get("_id")));
        System.out.println(feature.toString());
        System.out.println("Id is: " + feature.getUUID());
        testObjBase(feature);
    }

    private static BSONObject createFeature() {
        BSONObject afeat = new BasicDBObject();
            afeat.put("name", "mRFP1");
            afeat.put("sequence", "atgcatcatcatcatcatcatta");
        
        return afeat;
    }
   
    
    private static Schema createSchema(){
        
        Set<ClothoField> thefields = newHashSet(new ClothoField("name", String.class, "GFPmut3", "the feature name", null, false, Access.PUBLIC),
                                                new ClothoField("sequence", String.class, "atgcatgagatcatgcagccaactatttattaa", "the feature sequence", null, false, Access.PUBLIC));

        Schema schema = new ClothoSchema("Feature", "Act ontology standard representation ofa genetic feature", null, null, thefields);
        System.out.println(schema.toString());
        return schema;
    }
    
    private static boolean testSchema(Schema schema) {
        try {
            //Convert to JSON
//            JSONObject json = schema.toJSON();
            //Convert back to Schema
            JavaSchema converted = null; //TODO: JavaSchema.deserialize(json.toString());

//            System.out.println(converted.toJSON()); 
            
            
            if(!converted.getUUID().equals(schema.getUUID())) {
                return false;
            }
            
            
            System.out.println("Schema looks to be fine!!!");
            return true;
        } catch(Exception err) {
            err.printStackTrace();
            return false;  
        }
    }

    private static boolean testObjBase(ObjBase obj) {
        try {
            //Convert to JSON
//            JSONObject json = obj.toJSON();
            //Convert back to Schema
            ObjBase converted = null; //TODO: deserialize here 
            
            
            if(obj.getUUID().equals(converted.getUUID())) {
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

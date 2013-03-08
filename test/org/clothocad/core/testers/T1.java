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

import java.util.ArrayList;
import java.util.List;

import org.clothocad.core.aspects.Persistor;
import org.clothocad.core.datums.Instance;
import org.clothocad.core.datums.Schema;
import org.clothocad.core.datums.objbases.Person;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.datums.util.FieldType;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.datums.util.ServerScript;
import org.clothocad.core.util.Logger;
import org.json.JSONObject;

/**
 * Create a schema to represent a DNA sequence, and then
 * create a GFPuv instance of that SimpleFeature schema.
 * 
 * @author John Christopher Anderson
 */
public class T1 {
    public static void main(String[] args) {
        Schema feature = createFeatureSchema();
        Instance gfp = createGFPInstance(feature);
        
        if(feature==null || gfp == null)  {
            System.out.println("!!!!!!!!!!!!!!  Something wrong in T1 !!!!!!!!!!!!!!!  ");
            return;
        }
        System.out.println("T1:  All is good!");
    }
    
    private static Schema createFeatureSchema() {
        try {
            //Construct the "EmailAddress" schema and it's String content
            ServerScript indexer = new ServerScript("return true;", Language.JavaScript);
            ServerScript queryer = new ServerScript("return true;", Language.JavaScript);
            ServerScript validator = new ServerScript(""
                    //+ "clotho.say('+++person validator starting:   ' + inputs.email);"
                    + "clotho.log(\"INFO\", '+++feature validator starting:   ' + inputs.sequence);"
                    + "return true;"
                    + "", Language.JavaScript);

            List<ClothoField> thefields = new ArrayList<ClothoField>();

            ClothoField name = new ClothoField("name", FieldType.NAME, "mRFP13", 0);
            thefields.add(name);
            
            ClothoField seq = new ClothoField("sequence", FieldType.STRING, "ATGAGGAGAGGCATAGATTCAGCACCATGACCACCATGCAGAGTAA", 0);
            thefields.add(seq);

            Schema featureSchema = Schema.create(
                            Person.getAdmin(), 
                            "SimpleFeature", 
                            "A simple and sloppy representation of a Feature or other DNA sequence",
                            thefields, indexer, queryer, validator);

            //Change the UUID
            JSONObject obj = featureSchema.toJSON();
            obj.put("id", "specific-simplefeature-is-uuid");

            featureSchema = Schema.deserialize(obj.toString());
            
            // Chris' version:
            Persistor.get().persist(featureSchema);
            
            // Stephanie's version:
            //featureSchema.save();
            
            Logger.log(Logger.Level.INFO, featureSchema.getId() + "  " + featureSchema.getName() + "\n" + featureSchema.getDescription() + "\n...was created successfully, all good!");
            return featureSchema;
        } catch (Exception ex) {
            ex.printStackTrace();
            return null;
        }
    }
    
    private static Instance createGFPInstance(Schema feature) {
        try {
            //Create the Instance
            JSONObject fields = new JSONObject();
            
            fields.put("name",  "GFPuv" );
            fields.put("sequence",  "ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA" );

            String str = fields.toString();
            Instance gfp = Instance.create(Person.getAdmin(), feature, str);

            //Change the Id
            JSONObject obj = gfp.toJSON();
            obj.put("id", "specific-gfpuv-is-uuid");
            gfp = Instance.deserialize(obj.toString());
            
            // Chris' version:
            Persistor.get().persist(gfp);
            
            // Stephanie's version:
            //gfp.save();
            
            Logger.log(Logger.Level.INFO, gfp.getId() + "   " + gfp.getString("name") + "\n" + gfp.getString("sequence") + "\n...was created successfully, all good!");
            return gfp;
        } catch (Exception ex) {
            ex.printStackTrace();
            return null;
        }
    }
    
}

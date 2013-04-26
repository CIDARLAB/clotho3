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

import com.github.jmkgreen.morphia.logging.MorphiaLoggerFactory;
import com.github.jmkgreen.morphia.logging.slf4j.SLF4JLogrImplFactory;
import com.google.common.collect.Sets;
import com.mongodb.BasicDBObject;
import java.net.UnknownHostException;
import java.util.Set;
import org.bson.BSONObject;
import org.bson.types.ObjectId;
import org.clothocad.core.aspects.Persistor;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.layers.persistence.DBClassLoader;
import org.clothocad.core.layers.persistence.mongodb.MongoDBConnection;
import org.clothocad.core.schema.ClothoSchema;
import org.clothocad.core.schema.Schema;
import org.junit.Test;
import static org.objectweb.asm.Opcodes.*;

import static org.junit.Assert.*;
import org.junit.BeforeClass;

/**
 * Create a schema to represent a DNA sequence, and then
 * create a GFPuv instance of that SimpleFeature schema.
 * 
 * @author John Christopher Anderson
 */
public class Schema1Test {
    
    
    static {
    MorphiaLoggerFactory.registerLogger(SLF4JLogrImplFactory.class);
    }
    
    @BeforeClass
    public static void setUpClass() throws UnknownHostException {
        p.connect();
    }
    
    static Persistor p = new Persistor(new MongoDBConnection());
    static DBClassLoader cl = new DBClassLoader(p);
    
    private Schema createFeatureSchema() {
            Set<ClothoField> fields = Sets.newHashSet(new ClothoField("sequence", String.class, "ATACCGGA", "the sequence of the feature", null, false, ACC_PUBLIC));
        
            ClothoSchema featureSchema = new ClothoSchema("SimpleFeature", "A simple and sloppy representation of a Feature or other DNA sequence", null, null, fields);

            ObjectId id = new ObjectId();
            featureSchema.setUUID(id);
            p.save(featureSchema);
            
            return p.get(ClothoSchema.class, id);
    }
    
    
    @Test
    public void testClothoSchemaInstantiate() throws ClassNotFoundException, NoSuchFieldException, IllegalArgumentException, IllegalAccessException {
            Schema featureSchema = createFeatureSchema();
        
            BSONObject data = new BasicDBObject();
            
            String sequence = "ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA";
            
            data.put("name",  "GFPuv" );
            data.put("sequence",  sequence) ; //"ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA" );
            ObjectId id = new ObjectId();
            data.put("_id", id);
            
            p.save(data);
            
            ObjBase featureInstance = p.get(featureSchema.getEnclosedClass(cl), id);
            
            assertEquals("GFPuv", featureInstance.getName());
            assertEquals(sequence, featureInstance.getClass().getDeclaredField("sequence").get(featureInstance));
    }
    
    @Test
    public void testClothoSchemaCompile() throws ClassNotFoundException, NoSuchFieldException {
        Schema featureSchema = createFeatureSchema();
        Class featureClass = cl.loadClass(featureSchema.getBinaryName());
        
        assertEquals(ObjBase.class, featureClass.getSuperclass());
        
        assertEquals(2, featureClass.getDeclaredFields().length);
        //SCHEMA_NAME and sequence are the declared fields
        assertNotNull(featureClass.getDeclaredField("sequence"));
        
    }
    
}

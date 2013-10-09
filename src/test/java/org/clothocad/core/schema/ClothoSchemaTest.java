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
package org.clothocad.core.schema;

import com.google.common.collect.Sets;
import com.google.inject.Injector;
import com.mongodb.BasicDBObject;
import java.net.UnknownHostException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import javax.validation.ConstraintViolation;
import javax.validation.ConstraintViolationException;
import javax.validation.Validation;
import javax.validation.Validator;
import javax.validation.constraints.Pattern;
import org.bson.BSONObject;
import org.bson.types.ObjectId;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.persistence.DBClassLoader;
import org.clothocad.core.schema.Access;
import org.clothocad.core.schema.ClothoSchema;
import org.clothocad.core.schema.Constraint;
import org.clothocad.core.schema.Schema;
import org.clothocad.core.util.TestUtils;
import org.junit.Test;

import static org.junit.Assert.*;
import org.junit.BeforeClass;

public class ClothoSchemaTest {
    //TODO: 
    // validation
    // reference to other class

    static {
        //MorphiaLoggerFactory.registerLogger(SLF4JLogrImplFactory.class);
    }

    @BeforeClass
    public static void setUpClass() throws UnknownHostException {

        Injector injector = TestUtils.getDefaultTestInjector();
        p = injector.getInstance(Persistor.class);
        cl = injector.getInstance(DBClassLoader.class);

        p.connect();
        p.deleteAll();
        featureSchema = createFeatureSchema();

    }
    static Persistor p;
    static DBClassLoader cl;
    static Validator validator = Validation.buildDefaultValidatorFactory().getValidator();
    static Schema featureSchema;

    public static Schema createFeatureSchema() {

        ClothoField field = new ClothoField("sequence", String.class, "ATACCGGA", "the sequence of the feature", false, Access.PUBLIC);
        field.setConstraints(Sets.newHashSet(new Constraint("pattern", "regexp", "[ATUCGRYKMSWBDHVN]*", "flags", new Pattern.Flag[]{Pattern.Flag.CASE_INSENSITIVE})));
        Set<ClothoField> fields = Sets.newHashSet(field);

        ClothoSchema featureSchema = new ClothoSchema("SimpleFeature", "A simple and sloppy representation of a Feature or other DNA sequence", null, null, fields);

        ObjectId id = new ObjectId();
        featureSchema.setUUID(id);
        p.save(featureSchema);

        return p.get(ClothoSchema.class, id);
    }

    private ObjBase instantiateSchema(BSONObject data, Schema schema) throws ClassNotFoundException {
        ObjectId id = new ObjectId();
        data.put("_id", id);

        p.save(data.toMap());

        return p.get(schema.getEnclosedClass(cl), id);
    }

    @Test
    public void testClothoSchemaInstantiate() throws ClassNotFoundException, NoSuchFieldException, IllegalArgumentException, IllegalAccessException {

        BSONObject data = new BasicDBObject();

        String sequence = "ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA";

        data.put("name", "GFPuv");
        data.put("sequence", sequence); //"ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA" );

        ObjBase featureInstance = instantiateSchema(data, featureSchema);

        assertEquals("GFPuv", featureInstance.getName());
        assertEquals(sequence, featureInstance.getClass().getDeclaredField("sequence").get(featureInstance));
    }

    @Test
    public void testClothoSchemaValidate() throws ClassNotFoundException {
        //persistor.get auto-validates
        
        BSONObject data = new BasicDBObject();
        data.put("name", "BadSequence");
        data.put("sequence", "This is not a valid sequence.");
        try {
            ObjBase featureInstance = instantiateSchema(data, featureSchema);
            fail();
        } catch (ConstraintViolationException e) {

            Set<ConstraintViolation<?>> cvs = e.getConstraintViolations();
            assertTrue(cvs.size() == 1);
            ConstraintViolation<?> violation = cvs.iterator().next();
            assertEquals(violation.getMessage(), "must match \"[ATUCGRYKMSWBDHVN]*\"");
        }


        String sequence = "atgaGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATT"
                + "AGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACA"
                + "TACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAA"
                + "CACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAA"
                + "ACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCT"
                + "TTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTG"
                + "TTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAA"
                + "ACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGA"
                + "ATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACC"
                + "ATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCT"
                + "GTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTT"
                + "GAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA";

        data.put("name", "GFPuv");
        data.put("sequence", sequence);
        ObjBase featureInstance = instantiateSchema(data, featureSchema);
    }

    @Test
    public void testClothoSchemaCompile() throws ClassNotFoundException, NoSuchFieldException {

        Class featureClass = cl.loadClass(featureSchema.getBinaryName());

        assertEquals(ObjBase.class, featureClass.getSuperclass());

        assertEquals(2, featureClass.getDeclaredFields().length);
        //SCHEMA_NAME and sequence are the declared fields
        assertNotNull(featureClass.getDeclaredField("sequence"));

    }

    @Test
    public void testSchemaJSON() {
        Map output = p.toJSON(featureSchema);

        assertFalse(output.containsKey("isDeleted"));
        assertFalse(output.containsKey("lastUpdated") || output.containsKey("lastAccessed"));
        assertTrue(output.containsKey("schema"));
        assertEquals(Language.JSONSCHEMA.name(), output.get("language"));
        List fields = ((List) output.get("fields"));
        assertEquals(1, fields.size());
        Map field = (Map) fields.get(0);
        assertEquals("sequence", field.get("name"));
        assertNotNull(((Map) field.get("constraints")).get("pattern"));

        p.delete(featureSchema.getUUID());

        ObjectId id = p.save(p.toJSON(featureSchema));
        Schema secondSchema = p.get(ClothoSchema.class, featureSchema.getUUID());

        assertEquals(output, p.toJSON(secondSchema));
    }

    @Test
    public void testInstanceToJSON() throws ClassNotFoundException {
        BSONObject data = new BasicDBObject();

        String sequence = "atgaGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATT"
                + "AGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACA"
                + "TACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAA"
                + "CACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAA"
                + "ACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCT"
                + "TTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTG"
                + "TTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAA"
                + "ACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGA"
                + "ATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACC"
                + "ATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCT"
                + "GTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTT"
                + "GAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA";

        data.put("name", "GFPuv");
        data.put("sequence", sequence);

        ObjBase featureInstance = instantiateSchema(data, featureSchema);
        Map output = p.toJSON(featureInstance);

        assertNotNull(output);
        assertEquals(sequence, output.get("sequence"));
        assertEquals(data.get("name"), output.get("name"));
    }
}

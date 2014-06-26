/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.clothocad.model;

import com.fasterxml.jackson.core.JsonParseException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.schema.Converter;
import static org.clothocad.core.schema.ConverterTest.basicPartSchema;
import static org.clothocad.core.schema.ConverterTest.p;
import org.clothocad.core.schema.InferredSchema;
import org.clothocad.core.schema.Schema;
import org.clothocad.core.util.JSON;
import org.clothocad.core.util.TestUtils;
import org.junit.After;
import org.junit.AfterClass;
import static org.junit.Assert.*;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 *
 * @author prashantvaidyanathan
 */
public class SimpleSequenceTest {
    
    public SimpleSequenceTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    /**
     * Test of fromNucSeq method, of class SimpleSequence.
     */
  
    /**
     * Test of getSequence method, of class SimpleSequence.
     */
    
    
    
    
    
    
    
    
    /*
    public static final Persistor p = new TestUtils().getA(Persistor.class);
    public static Schema basicPartSchema = p.get(Schema.class, new ObjectId("org.clothocad.model.BasicPart"));

    public void testCanConvert() {
        Converter converter = new BasicPartConverter(p);
        Schema eugenePartSchema = new InferredSchema("eugene.dom.components.Part");
        assertTrue(converter.canConvert(eugenePartSchema));
    }

    @Test
    public void testConvertsTo() {
        Converter converter = new BasicPartConverter(p);
        assertEquals(basicPartSchema, converter.convertsTo());
    }

    @Test
    public void testConvert() throws JsonParseException, IOException {
        Converter<BasicPart> converter = new BasicPartConverter(p);
        Schema eugenePartSchema = new InferredSchema("eugene.dom.components.Part");
        Map<String, Object> eugeneJSON = JSON.deserializeObjectToMap("    {\n"
                + "         \"Name\":\"B0015\",\n"
                + "         \"schema\":\"eugene.dom.components.Part\",\n"
                + "         \"PartType\":\"Terminator\",\n"
                + "         \"Sequence\":\"CCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGGTGAACGCTCTCTACTAGAGTCACACTGGCTCACCTTCGGGTGGGCCTTTCTGCGTTTATA\",\n"
                + "         \"Pigeon\":\"t B0015\"\n"
                + "      }");

        BasicPart convertedPart = converter.convert(eugeneJSON, eugenePartSchema);
        assertEquals(Part.PartFunction.TERMINATOR, convertedPart.getType());
        assertEquals(eugeneJSON.get("Sequence").toString(), convertedPart.getSequence().getSeq());
        assertEquals(eugeneJSON.get("Name").toString(), convertedPart.getName());
    }
    */
    
    public static final Persistor p = new TestUtils().getA(Persistor.class);
    public static Schema nseqSchema = p.get(Schema.class, new ObjectId("org.clothocad.model.NucSeq"));
    
    public void testCanConvert() 
    {
        Converter converter = new NucSeqConverter(p);
        Schema simpleSeqSchema = new InferredSchema("org.clothocad.model.SimpleSequence");
        if(converter.canConvert(simpleSeqSchema))
        {
            System.out.println("Yes. SimpleSeq can be converted to a NucSeq");
        }
        else
            System.out.println(" SimpleSeq Cannot be converted to a NucSeq!!");
            
    }


    public void testConvertsTo() {
        Converter converter = new NucSeqConverter(p);
        if(nseqSchema.equals(converter.convertsTo()))
            System.out.println("ConvertsTo worked!");
        else
            System.out.println("ConvertsTo Failed!!");    
        //assertEquals(nseqSchema, converter.convertsTo());
    }

   
    public void testConvert() throws JsonParseException, IOException 
    {
        Converter<NucSeq> converter = new NucSeqConverter(p);
        Schema simpleSeqSchema = new InferredSchema("org.clothocad.model.SimpleSequence");
        Map<String, Object> simlpleSjson = new HashMap<String,Object>();
        
        

        simlpleSjson.put("Name","Just_Another_SimpleSeq");
        simlpleSjson.put("schema","org.clothocad.model.SimpleSequence");
        simlpleSjson.put("sequence","CCAGGCATCARATAAA");
        
                /*JSON.deserializeObjectToMap("    {\n"
                + "         \"Name\":\"Its a Simple Sequence\",\n"
                + "         \"schema\":\"org.clothocad.model.SimpleSequence\",\n"
               + "          \"sequence\":\"CCAGGCATCAAATAAA\",\n"
                + "      }");*/
        
        System.out.println(simlpleSjson.size());
        NucSeq convertedPart = converter.convert(simlpleSjson, simpleSeqSchema);
        System.out.println("NucSeq Sequence:" + convertedPart.getSeq());
        if(convertedPart.isDegenerate())
        {
            System.out.println("Degenerate Sequence");
        }
        
        
        Map<String,Object> nucseqmap = new HashMap<String,Object>();
        try {
            nucseqmap =  convertedPart.getNucSeqMap();
        } catch (IllegalArgumentException ex) {
            Logger.getLogger(SimpleSequenceTest.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            Logger.getLogger(SimpleSequenceTest.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
        System.out.println("Iterator Map: ");
        
        Iterator nsit = nucseqmap.entrySet().iterator();
        while(nsit.hasNext())
        {
            Map.Entry pairs = (Map.Entry)nsit.next();
            System.out.println(pairs.getKey() +" : "+pairs.getValue());
        }
        
        System.out.println("sequence : " +nucseqmap.get("sequence").toString());
        
    }
    
    
    @Test
    public void testGetSequence() {
        System.out.println("getSequence");
        SimpleSequence instance = null;
        String expResult = "";
        String result = instance.getSequence();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of setSequence method, of class SimpleSequence.
     */
    @Test
    public void testSetSequence() {
        System.out.println("setSequence");
        String sequence = "";
        SimpleSequence instance = null;
        instance.setSequence(sequence);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }
    
}

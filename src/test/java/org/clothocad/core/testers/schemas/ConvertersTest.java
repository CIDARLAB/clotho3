/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.testers.schemas;

import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.schema.Converter;
import org.clothocad.core.schema.Converters;
import org.clothocad.core.schema.InferredSchema;
import org.clothocad.core.schema.Schema;
import org.clothocad.core.utils.TestUtils;
import org.clothocad.model.BasicPartConverter;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author spaige
 */
public class ConvertersTest {

    private static Persistor persistor = TestUtils.getA(Persistor.class); 

    public ConvertersTest() {
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
    // TODO add test methods here.
    // The methods must be annotated with annotation @Test. For example:
    //
    // @Test
    // public void hello() {}

    private Converters prepareConverters() {
        Converters converters = new Converters();
        converters.addConverter(new BasicPartConverter(persistor));
        return converters;
    }

    @Test
    public void testGetConverterSchemas() {
        Converters converters = prepareConverters();
        Schema basicPartSchema = persistor.get(Schema.class, persistor.resolveSelector("BasicPart", false));
        Iterable<Schema> convertibleSchemas = converters.getConverterSchemas(basicPartSchema);
        assertEquals(new InferredSchema("eugene.dom.components.Part"), convertibleSchemas.iterator().next());
    }

    @Test
    public void testGetConverter() {
        Converters converters = prepareConverters();
        Schema basicPartSchema = persistor.get(Schema.class, persistor.resolveSelector("BasicPart", false));
        Schema eugeneSchema = new InferredSchema("org.eugene.dom.Part");
        Converter converter = converters.getConverter(basicPartSchema, eugeneSchema);
        assertTrue(converter.canConvert(eugeneSchema));
        assertEquals(basicPartSchema, converter.convertsTo());
    }
}
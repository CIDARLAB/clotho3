/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.schema;

import com.google.inject.Injector;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.schema.Converter;
import org.clothocad.core.schema.Converters;
import org.clothocad.core.schema.InferredSchema;
import org.clothocad.core.schema.Schema;
import org.clothocad.core.util.TestUtils;
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

    public static Persistor persistor;
    {
        Injector injector = TestUtils.getDefaultTestInjector();
        persistor = injector.getInstance(Persistor.class);
    } 

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
        Schema compositePartSchema = persistor.get(Schema.class, persistor.resolveSelector("CompositePart", false));
        Iterable<Schema> convertibleSchemas = converters.getConverterSchemas(basicPartSchema);
        assertEquals(new InferredSchema("eugene.dom.components.Part"), convertibleSchemas.iterator().next());
        
        convertibleSchemas = converters.getConverterSchemas(compositePartSchema);
        assertFalse(convertibleSchemas.iterator().hasNext());
    }

    @Test
    public void testGetConverter() {
        Converters converters = prepareConverters();
        Schema basicPartSchema = persistor.get(Schema.class, persistor.resolveSelector("BasicPart", false));
        Schema eugeneSchema = new InferredSchema("eugene.dom.components.Part");
        Converter converter = converters.getConverter(eugeneSchema, basicPartSchema);
        assertTrue(converter.canConvert(eugeneSchema));
        assertEquals(basicPartSchema, converter.convertsTo());
    }
}
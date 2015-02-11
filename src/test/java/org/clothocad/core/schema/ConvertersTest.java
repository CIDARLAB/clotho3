/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.schema;

import com.google.inject.Injector;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.util.AuthorizedShiroTest;
import org.clothocad.core.util.TestUtils;
import org.clothocad.model.BasicPartConverter;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author spaige
 */
public class ConvertersTest extends AuthorizedShiroTest {

    public static Persistor persistor;
    {
        Injector injector = TestUtils.getDefaultTestInjector();
        persistor = injector.getInstance(Persistor.class);
        persistor.initializeBuiltInSchemas();
    } 

    public ConvertersTest() {
    }

    private Converters prepareConverters() {
        Converters converters = new Converters();
        converters.addConverter(new BasicPartConverter(persistor));
        return converters;
    }

    @Test
    public void testGetConverterSchemas() {
        Converters converters = prepareConverters();
        Schema partSchema = persistor.get(Schema.class, new ObjectId("org.clothocad.model.Part"));
//        Schema compositePartSchema = persistor.get(Schema.class, new ObjectId("org.clothocad.model.CompositePart"));
        Iterable<Schema> convertibleSchemas = converters.getConverterSchemas(partSchema);
        assertEquals(new InferredSchema("eugene.dom.components.Part"), convertibleSchemas.iterator().next());
        
//        convertibleSchemas = converters.getConverterSchemas(compositePartSchema);
//        assertFalse(convertibleSchemas.iterator().hasNext());
    }

    @Test
    public void testGetConverter() {
        Converters converters = prepareConverters();
        Schema partSchema = persistor.get(Schema.class, new ObjectId("org.clothocad.model.Part"));
        Schema eugeneSchema = new InferredSchema("eugene.dom.components.Part");
        Converter converter = converters.getConverter(eugeneSchema, partSchema);
        assertTrue(converter.canConvert(eugeneSchema));
        assertEquals(partSchema, converter.convertsTo());
    }
}
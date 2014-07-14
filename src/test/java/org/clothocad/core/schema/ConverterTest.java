/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.schema;

import com.fasterxml.jackson.core.JsonParseException;
import java.io.IOException;
import java.util.Map;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.util.AuthorizedShiroTest;
import org.clothocad.core.util.JSON;
import org.clothocad.model.BasicPart;
import org.clothocad.model.BasicPartConverter;
import org.clothocad.model.Part.PartFunction;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author spaige
 */
public class ConverterTest extends AuthorizedShiroTest {

    public ConverterTest() {
        super();
        p = injector.getInstance(Persistor.class);
        basicPartSchema = p.get(Schema.class, new ObjectId("org.clothocad.model.BasicPart"));
    }
    public Persistor p;
    public Schema basicPartSchema ;

    @Test
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
        assertEquals(PartFunction.TERMINATOR, convertedPart.getType());
        assertEquals(eugeneJSON.get("Sequence").toString(), convertedPart.getSequence().getSeq());
        assertEquals(eugeneJSON.get("Name").toString(), convertedPart.getName());
    }
}
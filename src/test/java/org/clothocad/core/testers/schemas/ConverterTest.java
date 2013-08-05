/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.testers.schemas;

import com.fasterxml.jackson.core.JsonParseException;
import java.util.Map;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.schema.Converter;
import org.clothocad.core.schema.InferredSchema;
import org.clothocad.core.schema.Schema;
import org.clothocad.core.util.JSON;
import org.clothocad.core.utils.TestUtils;
import org.clothocad.model.BasicPart;
import org.clothocad.model.BasicPartConverter;
import org.clothocad.model.Part.PartFunction;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author spaige
 */
public class ConverterTest {

    public ConverterTest() {
    }
    public static final Persistor p = new TestUtils().getA(Persistor.class);
    public static Schema basicPartSchema = p.get(Schema.class, p.resolveSelector("BasicPart", false));

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
    public void testConvert() throws JsonParseException {
        Converter<BasicPart> converter = new BasicPartConverter(p);
        Schema eugenePartSchema = new InferredSchema("eugene.dom.components.Part");
        Map<String, Object> eugeneJSON = JSON.deserializeObject("    {\n"
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
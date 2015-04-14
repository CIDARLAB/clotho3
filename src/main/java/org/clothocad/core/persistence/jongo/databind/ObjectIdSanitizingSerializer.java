/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jongo.databind;

import com.fasterxml.jackson.core.JsonGenerationException;
import com.fasterxml.jackson.core.JsonGenerator;
import com.fasterxml.jackson.databind.SerializerProvider;
import com.fasterxml.jackson.databind.ser.std.StdSerializer;
import java.io.IOException;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.jongo.MongoUtil;

/**
 *
 * @author spaige
 */
public class ObjectIdSanitizingSerializer extends StdSerializer<ObjectId> {
    public ObjectIdSanitizingSerializer(){
        super(ObjectId.class);
    }
    
    @Override
    public void serialize(ObjectId value, JsonGenerator jgen, SerializerProvider provider) throws IOException, JsonGenerationException {
        String valueString = MongoUtil.sanitizeFieldName(value.toString());
        jgen.writeFieldName(valueString);
    }
    
}

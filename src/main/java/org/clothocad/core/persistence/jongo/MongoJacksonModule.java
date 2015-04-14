package org.clothocad.core.persistence.jongo;

import org.clothocad.core.persistence.jongo.databind.StringDesanitizingDeserializer;
import org.clothocad.core.persistence.jongo.databind.StringSanitizingSerializer;
import com.fasterxml.jackson.core.Version;
import com.fasterxml.jackson.databind.module.SimpleKeyDeserializers;
import com.fasterxml.jackson.databind.module.SimpleModule;
import com.fasterxml.jackson.databind.module.SimpleSerializers;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.jongo.databind.ObjectIdDesanitizingDeserializer;
import org.clothocad.core.persistence.jongo.databind.ObjectIdSanitizingSerializer;

/**
 * Holds mongodb-specific jackson configuration
 * 
 * escapes dots in field names to full-width-period character
 * 
 * having one serializer/deserializer pair for each key that might produce a 
 * dotted field name is not ideal, but it is a workaround until jackson
 * implements context-sensitive character escapes
 * 
 * @author spaige
 */
public class MongoJacksonModule extends SimpleModule {

    public MongoJacksonModule() {
        super("MongoModule", new Version(0, 0, 1, null, "org.clothocad", "clotho"));
    }

    @Override
    public void setupModule(SetupContext context) {
        SimpleSerializers keySerializers = new SimpleSerializers();
        keySerializers.addSerializer(String.class, new StringSanitizingSerializer());
        keySerializers.addSerializer(ObjectId.class, new ObjectIdSanitizingSerializer());
        context.addKeySerializers(keySerializers);
        
        SimpleKeyDeserializers keyDeserializers = new SimpleKeyDeserializers();
        keyDeserializers.addDeserializer(String.class, new StringDesanitizingDeserializer());
        keyDeserializers.addDeserializer(ObjectId.class, new ObjectIdDesanitizingDeserializer());
        context.addKeyDeserializers(keyDeserializers);
    }
}

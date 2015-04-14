/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jongo.databind;

import com.fasterxml.jackson.databind.DeserializationContext;
import com.fasterxml.jackson.databind.deser.std.StdKeyDeserializer;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.jongo.MongoUtil;

/**
 *
 * @author spaige
 */
public class ObjectIdDesanitizingDeserializer extends StdKeyDeserializer {

    public ObjectIdDesanitizingDeserializer(){
        super(ObjectId.class);
    }
    
    @Override
    protected ObjectId _parse(String key, DeserializationContext ctxt) throws Exception {
        return new ObjectId(MongoUtil.desanitizeFieldName(key));
    }
    
}

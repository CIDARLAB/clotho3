/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jongo.databind;

import com.fasterxml.jackson.databind.DeserializationContext;
import com.fasterxml.jackson.databind.deser.std.StdKeyDeserializer;
import org.clothocad.core.persistence.jongo.MongoUtil;

/**
 *
 * @author spaige
 */
public class StringDesanitizingDeserializer extends StdKeyDeserializer{

    public StringDesanitizingDeserializer(){
        super(String.class);
    }
    
    @Override
    protected String _parse(String key, DeserializationContext ctxt) throws Exception {
        return MongoUtil.desanitizeFieldName(key);
    }
    
}

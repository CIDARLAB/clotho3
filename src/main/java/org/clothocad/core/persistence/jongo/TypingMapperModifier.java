package org.clothocad.core.persistence.jongo;

import com.fasterxml.jackson.databind.ObjectMapper;
import org.jongo.marshall.jackson.configuration.MapperModifier;

/**
 *
 * @author spaige
 */
public class TypingMapperModifier implements MapperModifier {

    @Override
    public void modify(ObjectMapper mapper) {
        /*if you change the default typing, change the projection/query 
        * templates in JongoConnection.getGroupPermissions and 
        * JongoConnection.getUserPermissions*/
        mapper.enableDefaultTypingAsProperty(ObjectMapper.DefaultTyping.OBJECT_AND_NON_CONCRETE, "@class");
    }
    
}

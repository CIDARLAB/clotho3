/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import com.fasterxml.jackson.annotation.JsonAutoDetect;
import com.fasterxml.jackson.annotation.JsonAutoDetect.Visibility;
import com.fasterxml.jackson.core.JsonGenerationException;
import com.fasterxml.jackson.core.JsonParseException;
import com.fasterxml.jackson.core.Version;
import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.JsonMappingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.SerializationFeature;
import com.fasterxml.jackson.databind.module.SimpleModule;
import java.io.IOException;
import java.io.StringWriter;
import java.util.List;
import java.util.Map;

/**
 *
 * @author spaige
 */
public class JSON {
    
    public final static ObjectMapper mapper = new ObjectMapper();
    static {
        mapper.registerModule(new ClothoJacksonModule());
    }
    private final static TypeReference<Map<String, Object>>  stringToObject = new TypeReference<Map<String, Object>>(){};
    
    public static String serializeJSONMap(Map object){
        return serializeJSONMap(object, false);
    }
    
    public static String serializeJSONMap(Map object, boolean pretty){
        StringWriter writer = new StringWriter();
        if (pretty) mapper.enable(SerializationFeature.INDENT_OUTPUT);
        try {
            mapper.writeValue(writer, object);
                            return writer.toString();
        }  catch (JsonGenerationException | JsonMappingException ex) {
            throw new RuntimeException(ex);
        } catch (IOException ex) {
        }
        if (pretty) mapper.disable(SerializationFeature.INDENT_OUTPUT);
        return null;
    }
    
    public static String serialize(Object o){
        StringWriter writer = new StringWriter();
        try {
            mapper.writeValue(writer, o);
            return writer.toString();
        } catch (IOException ex) {
            throw new RuntimeException(ex);
        }
    }
    
    public static Map<String, Object> mappify(Object o){
        try {
                //XXX: ugh
            return deserializeObject(serialize(o));
        } catch (JsonParseException ex) {
            throw new RuntimeException(ex);
        }
    }
    
    public static Map<String, Object> deserializeObject(String json) throws JsonParseException {
        try {
            Map<String,Object> object = mapper.readValue(json, stringToObject);
                    return object;
        } catch (JsonMappingException ex) {
            throw new RuntimeException(ex);
        } catch (IOException ex) {
        } 
        return null;
    }
    
    public static List deserializeList(String json) throws JsonParseException {
        try {
            List object = mapper.readValue(json, List.class);
            return object;
        }  catch (JsonMappingException ex) {
            throw new RuntimeException(ex);
        }catch (IOException ex) {
        }
        
        return null;
    }
    
    
    static class ClothoJacksonModule extends SimpleModule {
        public ClothoJacksonModule(){
            //TODO: why is this deprecated?
            super("ClothoModule", new Version(0,0,1,null,"org.clothocad","clotho"));
            
        }

        @Override
        public void setupModule(SetupContext context) {
            context.setMixInAnnotations(Object.class, DisableGetters.class);
        }
        
    }
    
    @JsonAutoDetect(getterVisibility=Visibility.NONE)
    static abstract class DisableGetters {
        
    } 
    
}

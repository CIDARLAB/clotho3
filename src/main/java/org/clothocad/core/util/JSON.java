/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import com.fasterxml.jackson.annotation.JsonAutoDetect;
import com.fasterxml.jackson.annotation.JsonAutoDetect.Visibility;
import com.fasterxml.jackson.core.JsonGenerationException;
import com.fasterxml.jackson.core.JsonParseException;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.core.Version;
import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.JsonMappingException;
import com.fasterxml.jackson.databind.MappingIterator;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.ObjectReader;
import com.fasterxml.jackson.databind.SerializationFeature;
import com.fasterxml.jackson.databind.module.SimpleModule;
import java.io.File;
import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import javax.persistence.EntityNotFoundException;
import lombok.extern.slf4j.Slf4j;
import org.bson.types.ObjectId;
import org.clothocad.core.communication.ServerSideAPI;
import org.clothocad.core.persistence.Persistor;

/**
 *
 * @author spaige
 */
@Slf4j
public class JSON {
    
    public final static ObjectMapper mapper = new ObjectMapper();
    static {
        mapper.registerModule(new ClothoJacksonModule());
    }
    private final static TypeReference<Map<String, Object>>  stringToObject = new TypeReference<Map<String, Object>>(){};

    private static List<Map>
    parseJSONFiles(File[] paths) {
        ObjectReader reader = new ObjectMapper().reader(Map.class);
        List<Map> out = new ArrayList<Map>();
        for (File path : paths) {
            Iterator<Map> it;
            try {
                it = reader.readValues(path);
            } catch (JsonProcessingException ex) {
                log.warn("{}: not valid JSON", path.getAbsolutePath());
                continue;
            } catch (IOException ex) {
                log.warn("Could not open {}", path.getAbsolutePath());
                continue;
            }
            while (it.hasNext())
                out.add(it.next());
        }
        return out;
    }

    private static boolean
    persistorHasObject(Persistor persistor, Map object) {
        try {
            persistor.resolveSelector(object.get("name").toString(), false);
        } catch (EntityNotFoundException e) {
            return false;
        }
        return true;
    }

    private static boolean
    insertObjectIntoAPI(Map object, ServerSideAPI api, Persistor persistor,
                        boolean overwrite) {
        if (overwrite) {
            try {
                ObjectId result = api.create(object);
                if (result == null)
                    api.set(object);
            } catch (RuntimeException e) {
                return false;
            }
        } else {
            if (persistorHasObject(persistor, object))
                return true;
            try {
                ObjectId result = api.create(object);
            } catch (RuntimeException e) {
                return false;
            }
        }
        return true;
    }
    
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
        return serialize(o, false);
    }
    
    public static String serialize(Object o, boolean pretty){
        StringWriter writer = new StringWriter();
        if (pretty) mapper.enable(SerializationFeature.INDENT_OUTPUT);
        try {
            mapper.writeValue(writer, o);
            //if (pretty) mapper.disable(SerializationFeature.INDENT_OUTPUT);
            return writer.toString();
        } catch (IOException ex) {
           // if (pretty) mapper.disable(SerializationFeature.INDENT_OUTPUT);
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

    public static void
    importTestJSON(String path, Persistor persistor, boolean overwrite) {
        ServerSideAPI api = new DummyAPI(persistor);
        List<Map> objects = parseJSONFiles(new File(path).listFiles());
        for (Map obj : objects)
            if (!insertObjectIntoAPI(obj, api, persistor, overwrite))
                log.error("Could not load some files: {}", objects.toString());
    }
}

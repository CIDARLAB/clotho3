/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import com.fasterxml.jackson.annotation.JsonAutoDetect;
import com.fasterxml.jackson.annotation.JsonTypeInfo;
import com.fasterxml.jackson.core.JsonGenerationException;
import com.fasterxml.jackson.core.JsonParseException;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.core.Version;
import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.JsonMappingException;
import com.fasterxml.jackson.databind.MappingIterator;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.ObjectReader;
import com.fasterxml.jackson.databind.ObjectWriter;
import com.fasterxml.jackson.databind.SerializationFeature;
import static com.fasterxml.jackson.databind.SerializationFeature.FAIL_ON_EMPTY_BEANS;
import com.fasterxml.jackson.databind.annotation.JsonDeserialize;
import com.fasterxml.jackson.databind.module.SimpleModule;
import java.io.File;
import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import javax.persistence.EntityNotFoundException;
import javax.script.Bindings;
import javax.script.SimpleBindings;
import javax.validation.ConstraintViolation;
import javax.validation.Path.Node;
import javax.validation.metadata.ConstraintDescriptor;
import lombok.extern.slf4j.Slf4j;
import org.clothocad.core.aspects.Interpreter.RadixTrie.PatriciaTrie;
import org.clothocad.core.aspects.Interpreter.RadixTrie.Trie;
import org.clothocad.core.communication.ServerSideAPI;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.persistence.jackson.JSONViews;

/**
 *
 * @author spaige
 */
@Slf4j
public class JSON {

    //TODO: accept mapper configuration details from Guice, so sync'd w/ 
    //      persistor
    public final static ObjectMapper mapper = new ObjectMapper();

    static {
        mapper.registerModule(new ClothoJacksonModule());
        mapper.disable(FAIL_ON_EMPTY_BEANS);
        //mapper.disableDefaultTyping();
        //write types into serialized objects
        //mapper.enableDefaultTypingAsProperty(ObjectMapper.DefaultTyping.JAVA_LANG_OBJECT, "schema");
    }
    private final static TypeReference<Map<String, Object>> stringToObject = new TypeReference<Map<String, Object>>() {
    };

    public static String serializeJSONMapForExternal(Map object, boolean pretty) {
        StringWriter writer = new StringWriter();
        if (pretty) {
            mapper.enable(SerializationFeature.INDENT_OUTPUT);
        }
        try {
            mapper.writeValue(writer, object);
            return writer.toString();
        } catch (JsonGenerationException | JsonMappingException ex) {
            throw new RuntimeException(ex);
        } catch (IOException ex) {
        }
        if (pretty) {
            mapper.disable(SerializationFeature.INDENT_OUTPUT);
        }
        return null;
    }

    public static String serializeForExternal(Object o) {
        return serializeForExternal(o, false);
    }

    //Serialize w/ public view - omit db-only fields
    public static String serializeForExternal(Object o, boolean pretty) {
        ObjectWriter w = mapper.writerWithView(JSONViews.Public.class);
        if (pretty) {
            w = w.with(SerializationFeature.INDENT_OUTPUT);
        }
        try {
            return w.writeValueAsString(o);

        } catch (JsonProcessingException ex) {
            log.warn("Could not serialize object", ex);
            return ex.getMessage();
        }
    }

    public static String serialize(Object o) throws IOException {
        return serialize(o, false);
    }

    //Serialize w/ internal view - all fields included
    public static String serialize(Object o, boolean pretty) throws IOException {
        StringWriter writer = new StringWriter();
        if (pretty) {
            mapper.enable(SerializationFeature.INDENT_OUTPUT);
        }
        try {
            mapper.writeValue(writer, o);
            return writer.toString();
        } finally {
            //if (pretty) mapper.disable(SerializationFeature.INDENT_OUTPUT);
        }
    }

    public static Map<String, Object> mappify(Object o) {
        try {
            //XXX: ugh
            return deserializeObjectToMap(serialize(o));
        } catch (IOException ex) {
            throw new RuntimeException("Could not mappify "+o.toString()+": "+ex.getMessage());
        }
    }

    public static Map<String, Object> deserializeObjectToMap(String json) throws JsonParseException{
        try {
            Map<String, Object> object = mapper.readValue(json, stringToObject);
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
        } catch (JsonMappingException ex) {
            throw new RuntimeException(ex);
        } catch (IOException ex) {
        }

        return null;
    }

    public static class ClothoJacksonModule extends SimpleModule {

        public ClothoJacksonModule() {
            //TODO: why is this deprecated?
            super("ClothoModule", new Version(0, 0, 1, null, "org.clothocad", "clotho"));

        }

        @Override
        public void setupModule(SetupContext context) {
            context.setMixInAnnotations(Object.class, DisableGetters.class);
            context.setMixInAnnotations(Collection.class, DisableTypeInfo.class);
            context.setMixInAnnotations(Map.class, DisableTypeInfo.class);
//            context.setMixInAnnotations(Array.class, DisableTypeInfo.class);
            
            //Default types for interfaces unknown to Jackson
            context.setMixInAnnotations(Bindings.class, UseSimpleBindings.class);
            context.setMixInAnnotations(Trie.class, UsePatriciaTrie.class);
            
            //and it's safer to use public interfaces on some classes
            context.setMixInAnnotations(ConstraintViolation.class, UseDefaultAutoDetect.class);
            context.setMixInAnnotations(ConstraintDescriptor.class, UseDefaultAutoDetect.class);
            context.setMixInAnnotations(Node.class, UseDefaultAutoDetect.class);
        }
    }

    @JsonAutoDetect(fieldVisibility = com.fasterxml.jackson.annotation.JsonAutoDetect.Visibility.ANY,
            getterVisibility = com.fasterxml.jackson.annotation.JsonAutoDetect.Visibility.NONE,
            isGetterVisibility = com.fasterxml.jackson.annotation.JsonAutoDetect.Visibility.NONE,
            setterVisibility = JsonAutoDetect.Visibility.NONE)
    static abstract class DisableGetters {
    }

    @JsonTypeInfo(use = JsonTypeInfo.Id.NONE)
    static abstract class DisableTypeInfo {
    }
    
    @JsonDeserialize(as = SimpleBindings.class)
    static abstract class UseSimpleBindings{
        
    }
    
    @JsonDeserialize(as = PatriciaTrie.class)
    static abstract class UsePatriciaTrie {
        
    }

    @JsonAutoDetect(fieldVisibility = JsonAutoDetect.Visibility.DEFAULT, 
            getterVisibility = JsonAutoDetect.Visibility.DEFAULT,
            isGetterVisibility = JsonAutoDetect.Visibility.DEFAULT)
    static abstract class UseDefaultAutoDetect {};
    
    //XXX: move to test utils
    public static void importTestJSON(String path, Persistor persistor, boolean overwrite) {
        ServerSideAPI api = new DummyAPI(persistor);
        ObjectReader reader = new ObjectMapper().reader(Map.class);

        List<Map> objects = new ArrayList<>();
        for (File child : new File(path).listFiles()) {
            if (!child.getName().endsWith(".json")) {
                continue;
            }
            try {
                MappingIterator<Map> it = reader.readValues(child);
                while (it.hasNext()) {
                    objects.add(it.next());
                }

            } catch (JsonProcessingException ex) {
                log.warn("Could not process {} as JSON", child.getAbsolutePath());
            } catch (IOException ex) {
                log.warn("Could not open {}", child.getAbsolutePath());
            }
        }

        while (objects.size() > 0) {
            int prevSize = objects.size();
            List<Map> newObjects = new ArrayList();

            for (Map obj : objects) {
                if (!overwrite) {
                    try {
                        persistor.resolveSelector(obj.get("name").toString(), false);
                        continue;

                    } catch (EntityNotFoundException e) {
                    }
                }
                try {
                    ObjectId result = api.create(obj);

                    if (overwrite && result == null) {
                        api.set(obj);
                    }
                } catch (RuntimeException e) {
                    newObjects.add(obj);
                } catch (ClassCircularityError e) {
                    e.getMessage();
                }
            }

            objects = newObjects;
            if (objects.size() >= prevSize) {
                log.error("Could not load some files: {}", objects.toString());
                return;
            }

        }
    }
}

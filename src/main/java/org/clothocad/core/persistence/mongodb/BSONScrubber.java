/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.mongodb;

import com.github.jmkgreen.morphia.mapping.MappedClass;
import com.github.jmkgreen.morphia.mapping.MappedField;
import com.github.jmkgreen.morphia.mapping.Mapper;
import com.mongodb.DBRef;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.clothocad.core.persistence.DBOnly;
import org.clothocad.core.persistence.Replace;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 *
 * @author spaige
 */
public class BSONScrubber {

    public BSONScrubber(Mapper mapper) {
        this.mapper = mapper;
    }
    Mapper mapper;
    Logger logger = LoggerFactory.getLogger(BSONScrubber.class);
    public static final Set<String> internalFieldNames = new HashSet<>();

    {
        internalFieldNames.add("className");
        internalFieldNames.add("otherClasses");
        //XXX: demo hack - find way to designate added fields as DBOnly
        internalFieldNames.add(ClothoMappedField.VIRTUAL_PREFIX + "binaryName");
    }

    public Map<String, Object> scrub(Map<String, Object> input) {
        return scrub(input, null);
    }

    public Map<String, Object> scrub(Map<String, Object> input, Set<MappedClass> classHints) {
        Map<String, Object> output = new HashMap<>();
        Set<MappedClass> classSet = getClassSet(input);
        if (classHints != null) {
            classSet.addAll(classHints);
        }
        Set<String> excludes = getExcludedKeys(classSet);

        for (String key : input.keySet()) {
            String keyString;
            if (excludes.contains(key)) {
                continue;
            }
            if (key.startsWith(ClothoMappedField.VIRTUAL_PREFIX)) {
                keyString = key.substring(ClothoMappedField.VIRTUAL_PREFIX.length());
            } //rename _id to id
            else if (key.equals("_id")) {
                keyString = "id";
            } else {
                keyString = key;
            }
            if (isBasicValue(input.get(key))) //XXX
            {
                output.put(keyString, input.get(key));
                //transform references to just uuid
            } else if (extractReferenceId(input.get(key)) != null) {
                output.put(keyString, extractReferenceId(input.get(key)));
            } 
            //schema field
            else {
                output.put(keyString, scrub(input.get(key), getClassHints(keyString, classSet)));
            }
        }


        return output;
    }

    public List scrub(List input, Set<MappedClass> classHints) {
        List output = new ArrayList();
        for (Object item : input) {
            output.add(scrub(item, classHints));
        }

        return output;
    }

    public Object scrub(Object input, Set<MappedClass> classHints) {
        if (isBasicValue(input)) {
            return input;
        } else if (input instanceof Map) {
            return scrub((Map) input, classHints);
        } else if (input instanceof List) {
            return scrub((List) input, classHints);
        }
        logger.warn("Non-basic, non-Map, non-List object: {}", input);
        return input;
    }

    private Object extractReferenceId(Object reference) {
        if (reference instanceof Map && ((Map) reference).containsKey("$ref")) {
            return ((Map) reference).get("$id");
        }
        if (reference instanceof DBRef) {
            return ((DBRef) reference).getId();
        }
        return null;
    }

    public boolean isBasicValue(Object input) {
        return (input instanceof String || input instanceof Boolean || input == null || input instanceof Integer || input instanceof Number);
    }

    protected Set<MappedClass> getClassSet(Map<String, Object> input) {
        Set<MappedClass> output = new HashSet<>();
        if (input.containsKey("className")) {
            Map<String, MappedClass> mappedClasses = mapper.getMCMap();
            if (mappedClasses.containsKey((String) input.get("className"))) {
                output.add(mappedClasses.get((String) input.get("className")));
            }
        }
        if (input.containsKey("otherClasses")) {
            for (String className : forceCollection(input.get("otherClasses"))) {
            }
        }

        return output;
    }

    protected Set<String> getExcludedKeys(Set<MappedClass> classes) {
        Set<String> excludedFieldNames = new HashSet<>();
        excludedFieldNames.addAll(internalFieldNames);

        for (MappedClass c : classes) {
            addStoredFieldNames(c.getFieldsAnnotatedWith(DBOnly.class), excludedFieldNames);
            //addStoredFieldNames(c.getFieldsAnnotatedWith(Replace.class), excludedFieldNames);
        }

        return excludedFieldNames;
    }

    protected Set<MappedClass> getClassHints(String fieldName, Set<MappedClass> parentClasses) {
        Set<MappedClass> classHints = new HashSet<>();
        for (MappedClass c : parentClasses) {
            MappedField field = c.getMappedField(fieldName);
            if (field != null) {
                if (field.getSubClass() != null) {
                    classHints.add(mapper.getMappedClass(field.getSubClass()));
                } else {
                    classHints.add(mapper.getMappedClass(field.getType()));
                }
            }
        }
        return classHints;
    }

    protected static void addStoredFieldNames(Collection<MappedField> fields, Set<String> fieldNames) {
        for (MappedField f : fields) {
            fieldNames.add(f.getNameToStore());
        }
    }

    protected static Collection<String> forceCollection(Object input) {
        if (input instanceof Collection) {
            return (Collection) input;
        } else if (input instanceof String[]) {
            return Arrays.asList((String[]) input);
        } else {
            throw new IllegalArgumentException("Could not convert input to Collection");
        }
    }
}

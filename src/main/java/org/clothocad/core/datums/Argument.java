/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.datums;

import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.extern.slf4j.Slf4j;
import org.bson.types.ObjectId;
import org.clothocad.core.persistence.IdUtils;
import org.clothocad.core.schema.Schema;

/**
 *
 * @author spaige
 */
@NoArgsConstructor
@Slf4j
public class Argument {

    @Getter
    private String name;
    @Getter
    //XXX: @Replace(encoder = "jsonifyFieldType", decoder = "decodeFieldType")
    private Class type;

    public Argument(String name, Class type) {
        this.name = name;
        this.type = type;
    }

    public String jsonifyFieldType() {
        Class c = this.type;

        if (Schema.isSchemaClassName(c.getName())) ///XXX: fix for inner classes
        {
            return Schema.extractIdFromClassName(c.getName());
        }
        if (ObjBase.class.isAssignableFrom(c)) {
            return c.getSimpleName();
        }
        if (ObjectId.class.isAssignableFrom(c)) {
            return "id";
        }
        if (Date.class.isAssignableFrom(c)) {
            return "date";
        }
        if (String.class.isAssignableFrom(c) || c.equals(char.class)) {
            return "string";
        }
        if (Boolean.class.isAssignableFrom(c) || c.equals(boolean.class)) {
            return "boolean";
        }
        if (Number.class.isAssignableFrom(c)
                || c.equals(byte.class)
                || c.equals(short.class)
                || c.equals(int.class)
                || c.equals(long.class)
                || c.equals(float.class)
                || c.equals(double.class)) {
            return "number"; // String.format("number(%s)", c.getSimpleName());
        }
        if (c.isArray() || Collection.class.isAssignableFrom(c)) {
            //todo: parameterize array types;
            return "array";

        }
        if (Map.class.isAssignableFrom(c)) {
            return "object";
        }
        log.warn("Unable to jsonify field type {}", c.getName());
        return "object";
    }

    protected final static Map<String,Class> classMap;
    static {
        //todo: make immutable
        classMap = new HashMap<>();
        classMap.put("string", String.class);
        classMap.put("boolean", Boolean.class);
        classMap.put("number", Number.class);
        classMap.put("array", List.class);
        classMap.put("id", ObjectId.class);
        classMap.put("date", Date.class);
        classMap.put("object", Map.class);
    }
    
    public void decodeFieldType(Map object) {
        String s = (String) object.get("type");

        Class c = classMap.get(s.toLowerCase());
        if (c != null){
            type = c;
            return;
        }   
        ObjectId id;
        if (ObjectId.isValid(s)){
        id = new ObjectId(s);            
        } else {
            id = IdUtils.resolveSelector(s, true);
        }
        

        try {
            
            type = IdUtils.getClass(id);
        } catch (ClassNotFoundException ex) {
            throw new RuntimeException("Could not find schema: "+ s, ex);
        }
    }
}
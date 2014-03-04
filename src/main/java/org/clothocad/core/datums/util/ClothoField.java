/*
 * 
Copyright (c) 2010 The Regents of the University of California.
All rights reserved.
Permission is hereby granted, without written agreement and without
license or royalty fees, to use, copy, modify, and distribute this
software and its documentation for any purpose, provided that the above
copyright notice and the following two paragraphs appear in all copies
of this software.

IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.

THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
ENHANCEMENTS, OR MODIFICATIONS..
 */

package org.clothocad.core.datums.util;

import com.google.common.collect.ImmutableMap;
import java.lang.reflect.Field;
import java.lang.reflect.GenericArrayType;
import java.lang.reflect.ParameterizedType;
import java.lang.reflect.Type;
import java.lang.reflect.TypeVariable;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import lombok.Getter;
import lombok.Setter;
import org.clothocad.core.persistence.annotations.Reference;
import org.clothocad.core.persistence.annotations.ReferenceCollection;
import org.clothocad.core.persistence.annotations.Rename;
import org.clothocad.core.schema.Access;
import org.clothocad.core.schema.Constraint;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
/**
 * @author John Christopher Anderson
 */

@Getter
@Setter
public class ClothoField {
    
    private ClothoField() {}
    
    public ClothoField(Field field){
        
        if (field.getAnnotation(Rename.class) != null){
            name = field.getAnnotation(Rename.class).value();
        } else {
            name = field.getName();            
        }
        
        type = field.getType();
        reference = field.getAnnotation(Reference.class) != null;
        referenceCollection = field.getAnnotation(ReferenceCollection.class) != null;
        
        //TODO: access, validate
        
        //TODO: metadata
        //example
        //description
    }
    
    private static final Logger logger = LoggerFactory.getLogger(ClothoField.class);
    
    public ClothoField(String name, Class type, String example, String description, boolean reference, Access access) {
        this.name = name;
        this.type = type;
        this.example = example;
        this.access = access; 
        this.reference = reference;
        this.description = description;
    }
    

    private String name;
    
    private Class<?> type;
    private Type subtype;
    private String example;   //A string representation/explanation of an expected value
    private Access access;  
    private boolean reference;
    private boolean referenceCollection;
    
    private Set<Constraint>  constraints;
    
    public Map prettyPrintConstraints(){
        if (constraints == null) return null;
        Map<String, Map<String,Object>> output = new HashMap<>();
        for (Constraint constraint : constraints){
            Map<String, Object> constraintMap = new HashMap<>();
            for (String key : constraint.getValues()){
                constraintMap.put(key, constraint.getValue(key));
            }
            output.put(constraint.getConstraint(), constraintMap);
        }
        return output;
    }
    
    //metadata
    private String description;
    
    public String getSetterName(){
        return "set" + capitalize(name);
    }
    
    public String getGetterName(){
        if (this.type.equals(Boolean.class)) return "is" + capitalize(name);
        else return "get" + capitalize(name);
    }
    
    private static String capitalize(String s){
        if (s.length() == 0) return s;
        return s.substring(0,1).toUpperCase() + s.substring(1);
    }
    /*
    public String jsonifyFieldType(){
        Class c = this.type;
        
        if (Schema.isSchemaClassName(c.getName())) ///XXX: fix for inner classes
            return Schema.extractIdFromClassName(c.getName());
        if (ObjBase.class.isAssignableFrom(c)) return c.getSimpleName();
        if (ObjectId.class.isAssignableFrom(c)) return "id";
        if (Date.class.isAssignableFrom(c)) return "date";
        if (String.class.isAssignableFrom(c) || c.equals(char.class)) return "string";
        if (Boolean.class.isAssignableFrom(c) || c.equals(boolean.class)) return "boolean";
        if (Number.class.isAssignableFrom(c) || 
                c.equals(byte.class) || 
                c.equals(short.class) ||
                c.equals(int.class) ||
                c.equals(long.class) ||
                c.equals(float.class) ||
                c.equals(double.class)) return "number"; // String.format("number(%s)", c.getSimpleName());
        if (c.isArray() || Collection.class.isAssignableFrom(c)){
           //todo: parameterize array types;
            return "array";
            
        }
        if (Map.class.isAssignableFrom(c)) return "object";
        logger.warn("Unable to jsonify field type {}", c.getName());
        return "object";
    }
    
    public void decodeFieldType(Map object){
        String s = (String) object.get("javaType");
        if (s == null) s= (String) object.get(ClothoMappedField.VIRTUAL_PREFIX + "javaType");
        try {
            type = Class.forName(s, true, Schema.cl);
        } catch (ClassNotFoundException ex) {
            throw new RuntimeException(ex);
        }
    }
    */
    public void decodeConstraints(Map object){
        Map<String, Map<String, Object>> constraints = (Map<String, Map<String, Object>>) object.get("constraints");
        if (constraints == null) return;
        
        Set<Constraint> realConstraints = new HashSet<>();
        for (String constraint : constraints.keySet()){
            realConstraints.add(new Constraint(constraint, constraints.get(constraint)));
        }
        this.constraints = realConstraints;
    }
    
    public Class<?> getType(){
        return wrap(type);
    }
  
    //Morphia can't decode primitive classes
    //http://stackoverflow.com/questions/1704634/simple-way-to-get-wrapper-class-type-in-java
    // safe because both Long.class and long.class are of type Class<Long>
    @SuppressWarnings("unchecked")
    private static <T> Class<T> wrap(Class<T> c) {
        return c.isPrimitive() ? (Class<T>) PRIMITIVES_TO_WRAPPERS.get(c) : c;
    }
    private static final Map<Class<?>, Class<?>> PRIMITIVES_TO_WRAPPERS = new ImmutableMap.Builder<Class<?>, Class<?>>()
            .put(boolean.class, Boolean.class)
            .put(byte.class, Byte.class)
            .put(char.class, Character.class)
            .put(double.class, Double.class)
            .put(float.class, Float.class)
            .put(int.class, Integer.class)
            .put(long.class, Long.class)
            .put(short.class, Short.class)
            .put(void.class, Void.class)
            .build();
    //Constraints
    
    //#
    //multipleof
    //maximum
    //exclusivemaximum
    //minimum
    //exclusiveminimum
    
    //size
    //pattern (regex match)
    
    
    //notnull
    
    public static Type getParameterizedType(Type type) {
        int index = 0;
        if (type instanceof ParameterizedType) {
            ParameterizedType ptype = (ParameterizedType) type;
            if ((ptype.getActualTypeArguments() != null) && (ptype.getActualTypeArguments().length <= index)) {
                return null;
            }
            Type paramType = ptype.getActualTypeArguments()[index];
            if (paramType instanceof GenericArrayType) {
                return ((GenericArrayType) paramType).getGenericComponentType();
            } else {
                if (paramType instanceof ParameterizedType) {
                    return paramType;
                } else {
                    if (paramType instanceof TypeVariable) {
                        // TODO: Figure out what to do... Walk back up the to
                        // the parent class and try to get the variable type
                        // from the T/V/X
//						throw new MappingException("Generic Typed Class not supported:  <" + ((TypeVariable) paramType).getName() + "> = " + ((TypeVariable) paramType).getBounds()[0]);
                        return paramType;
                    } else if (paramType instanceof Class) {
                        return (Class) paramType;
                    } else {
                        throw new RuntimeException("Unknown type... pretty bad... call for help, wave your hands... yeah!");
                    }
                }
            }
        }
        return null;
    }
    
}
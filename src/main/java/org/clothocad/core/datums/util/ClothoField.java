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

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import lombok.Getter;
import lombok.Setter;
import org.clothocad.core.datums.Function;
import org.clothocad.core.layers.persistence.Add;
import org.clothocad.core.layers.persistence.Replace;
import org.clothocad.core.layers.persistence.mongodb.ClothoMappedField;
import org.clothocad.core.schema.Access;
import org.clothocad.core.schema.Constraint;
import org.clothocad.core.schema.Schema;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
/**
 * @author John Christopher Anderson
 */

@Getter
@Setter
@Add(provider="getType", name="javaType")
public class ClothoField {
    
    private ClothoField() {}
    
    private static final Logger logger = LoggerFactory.getLogger(ClothoField.class);
    
    public ClothoField(String name, Class type, String example, String description, Function validate, boolean reference, Access access) {
        this.name = name;
        this.type = type;
        this.example = example;
        this.access = access; 
        this.reference = reference;
        this.validate = validate;
        this.description = description;
    }
    

    private String name;
    
    @Replace(encoder="jsonifyFieldType", decoder="decodeFieldType")
    private Class<?> type;
    private String example;   //A string representation/explanation of an expected value
    private Access access;  
    private boolean reference;
    private Function validate;
    
    @Replace(encoder="prettyPrintConstraints", decoder="decodeConstraints")
    private Set<Constraint>  constraints;
    
    public Map prettyPrintConstraints(){
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
    
    public String jsonifyFieldType(){
        return jsonifyFieldType(type);
    }
    
    public static String jsonifyFieldType(Class c){
        if (Schema.isSchemaClassName(c.getName())) ///XXX: fix for inner classes
            return Schema.extractIdFromClassName(c.getName());
        if (String.class.isAssignableFrom(c)) return "string";
        if (Boolean.class.isAssignableFrom(c)) return "boolean";
        if (Number.class.isAssignableFrom(c)) return "number"; // String.format("number(%s)", c.getSimpleName());
        if (c.isArray() || Collection.class.isAssignableFrom(c)){
            return "array";
        }
        if (Map.class.isAssignableFrom(c)) return "object";
        logger.warn("Unable to jsonify field type {}", c.getName());
        return "?";
    }
    
    public void decodeFieldType(Map object){
        String s = (String) object.get("javaType");
        if (s == null) object.get(ClothoMappedField.VIRTUAL_PREFIX + "javaType");
        try {
            type = Class.forName(s, true, Schema.cl);
        } catch (ClassNotFoundException ex) {
            throw new RuntimeException(ex);
        }
    }
    
    public void decodeConstraints(Map object){
        Set<Constraint> realConstraints = new HashSet<>();
        Map<String, Map<String, Object>> constraints = (Map<String, Map<String, Object>>) object.get("constraints");
        for (String constraint : constraints.keySet()){
            realConstraints.add(new Constraint(constraint, constraints.get(constraint)));
        }
        this.constraints = realConstraints;
    }
    
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
    
    
    
}
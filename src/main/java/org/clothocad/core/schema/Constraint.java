/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.schema;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import javax.validation.constraints.AssertFalse;
import javax.validation.constraints.AssertTrue;
import javax.validation.constraints.Max;
import javax.validation.constraints.Min;
import javax.validation.constraints.NotNull;
import javax.validation.constraints.Pattern;
import javax.validation.constraints.Size;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.clothocad.core.schema.constraintutils.ConstraintValues;
import org.clothocad.core.schema.constraintutils.PatternValues;

/**
 *
 * @author spaige
 */
@NoArgsConstructor
public class Constraint {
    public Constraint(String constraint, Object... args){
        values = getConstraintValue(constraint, null);
        this.constraint = constraint;
        if (args.length % 2 != 0){
            throw new IllegalArgumentException("Needs value for each value name");
        }
        for (int i=0; i< args.length; i=i+2){
            if (!(args[i] instanceof String)){
                throw new IllegalArgumentException("Expected string value name");
            }
            values.put((String) args[i], args[i+1]);
        }
    }
    
    public Constraint(String constraint, Map<String, Object> values){
        this.constraint = constraint;
        this.values = getConstraintValue(constraint, values);
    }
    
    @Getter
    @Setter
    protected String constraint;
    
    private static final Map<String, Class> annotations;
    static {
        annotations = new HashMap<>();
        annotations.put("pattern", Pattern.class);
        annotations.put("true", AssertTrue.class);
        annotations.put("false", AssertFalse.class);
        annotations.put("max", Max.class);
        annotations.put("min", Min.class);
        annotations.put("notNull", NotNull.class);
        annotations.put("size", Size.class);
    }
    
    protected ConstraintValues values;
    
    public Set<String> getValues(){
       Set<String> set  =  values.keySet();
       return set;
    }
    
    public Object getValue(String key){
        return values.get(key);
    }
    
    public Class getAnnotation(){
        return annotations.get(constraint);
    }

    protected ConstraintValues getConstraintValue(String constraint, Map<String, Object> values){
        ConstraintValues cvs;
        switch(constraint){
            case("pattern"):
                cvs = new PatternValues();
                break;
            default:
                throw new IllegalArgumentException();
        }
        
        if (values != null){
            for (String key : values.keySet()){
                cvs.put(key, values.get(key));
            }
        }
        
        return cvs;
    }
}

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence;

import com.github.jmkgreen.morphia.annotations.Embedded;
import java.lang.annotation.Annotation;
import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

/**
 *
 * @author spaige
 */
@Target(value = {ElementType.FIELD})
@Retention(value = RetentionPolicy.RUNTIME)
public @interface Replace {
    public String value() default "";
    
    public String encoder() default "";
    public String decoder() default "";
    
    public Class concreteClass() default Object.class;
    
    public Class<? extends Annotation> type() default Embedded.class; 
    
    //doesn't do anything until per-field change tracking is implemented
    public String[] dependsOn() default {};
}

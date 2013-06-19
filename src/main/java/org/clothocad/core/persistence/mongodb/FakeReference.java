/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.mongodb;

import com.github.jmkgreen.morphia.annotations.Reference;
import com.github.jmkgreen.morphia.mapping.Mapper;
import java.lang.annotation.Annotation;

/**
 *
 * @author spaige
 */
public class FakeReference implements Annotation, Reference {
    //TODO: change from default values

    @Override
    public Class<? extends Annotation> annotationType() {
        return Reference.class;
    }

    @Override
    public String value() {
        return Mapper.IGNORED_FIELDNAME;
    }

    @Override
    public Class<?> concreteClass() {
        return Object.class;
    }

    @Override
    public boolean ignoreMissing() {
        return false;
    }

    @Override
    public boolean lazy() {
        return false;
    }
    
}

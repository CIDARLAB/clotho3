/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.mongodb;

import com.github.jmkgreen.morphia.mapping.MappedClass;
import com.github.jmkgreen.morphia.mapping.MappedField;
import java.lang.annotation.Annotation;
import java.util.ArrayList;
import java.util.List;
import org.clothocad.core.persistence.Add;
import org.clothocad.core.persistence.DBOnly;
import org.clothocad.core.persistence.Replace;

/**
 *
 * @author spaige
 */
class ClothoMappedClass extends MappedClass {

    public ClothoMappedClass(Class type, ClothoMapper aThis) {
        super(type, aThis);
    }

    @Override
    protected void discover() {
        super.discover();
        //make virtual fields
        List<Annotation> anns = this.getAnnotations(Add.class);
        if (anns != null) for (Annotation ann : anns){
            Add add = (Add) ann;
            this.getMappedFields().add(new ClothoMappedField(add, this.getClazz()));
        }
        
        List<MappedField> fields = this.getFieldsAnnotatedWith(Replace.class);
        if (fields != null) for (MappedField f : fields){
            this.getMappedFields().remove(f);
            this.getMappedFields().add(new ClothoMappedField(f.getAnnotation(Replace.class), f, this.getClazz()));
        }
    }
    
    
}

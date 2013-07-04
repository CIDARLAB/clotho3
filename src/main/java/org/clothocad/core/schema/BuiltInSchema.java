/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.schema;

import com.github.jmkgreen.morphia.utils.ReflectionUtils;
import java.lang.reflect.Field;
import java.util.HashSet;
import java.util.Set;
import lombok.NoArgsConstructor;
import lombok.extern.slf4j.Slf4j;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.model.Person;

/**
 *
 * @author spaige
 */
@NoArgsConstructor
@Slf4j
public class BuiltInSchema extends JavaSchema {

    public BuiltInSchema(Class<? extends ObjBase> c) {
        super();
        fields = extractFields(c);
        //TODO: functions
        
        binaryName = c.getCanonicalName();
        setName(c.getSimpleName());
        this.c = c;
    }

    Class<? extends ObjBase> c;
    
    @Override
    public void setSource(String source) {
        throw new UnsupportedOperationException("You cannot set the source of a built-in schema.");
    }

    private String binaryName;
    
    @Override
    public String getBinaryName() {
        return binaryName;
    }
    
    @Override
    public Class<? extends ObjBase> getEnclosedClass(ClassLoader cl){
        return c;
    }

    private Set<ClothoField> extractFields(Class c) {
        Set<ClothoField> fields = new HashSet<>();
        for (Field field : ReflectionUtils.getDeclaredAndInheritedFields(c, true)) {
            fields.add(new ClothoField(field));
        }
        return fields;
    }
}

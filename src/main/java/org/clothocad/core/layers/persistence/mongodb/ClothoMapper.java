/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.layers.persistence.mongodb;

import com.github.jmkgreen.morphia.mapping.DefaultMapper;
import com.github.jmkgreen.morphia.mapping.MappedClass;
import com.github.jmkgreen.morphia.mapping.MappedField;
import com.github.jmkgreen.morphia.mapping.MapperOptions;
import com.github.jmkgreen.morphia.mapping.lazy.proxy.ProxyHelper;
import com.mongodb.BasicDBObject;
import com.mongodb.DBObject;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import lombok.Getter;
import lombok.Setter;
import org.bson.BSONObject;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.layers.persistence.Add;
import org.clothocad.core.layers.persistence.DBOnly;
import org.clothocad.core.layers.persistence.Replace;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
/**
 *
 * @author spaige
 */

//TODO: perferentially remove schema classes from cache
import org.clothocad.core.aspects.JSONSerializer;
public class ClothoMapper extends DefaultMapper implements JSONSerializer{
    public ClothoMapper(MapperOptions opts){
        super(opts);
    }
    
    public ClothoMapper(){
        super(defaultOptions);
    }
    
    static final Logger logger = LoggerFactory.getLogger(ClothoMapper.class);
    static final MapperOptions defaultOptions = new MapperOptions();
    static {
        defaultOptions.objectFactory = new PolymorphicObjectFactory();
    }
    static {
        //add our annotations to the interesting annotations lists
        MappedClass.interestingAnnotations.add(Add.class);
        MappedField.interestingAnnotations.add(DBOnly.class);
        MappedField.interestingAnnotations.add(Replace.class);
    }
    @Getter  
    @Setter
    //XXX: circular reference (this -> cl -> persistor -> this)
    ClassLoader cl;
    

    
    public BSONObject toJSON(Map data){
        //TODO: cache
        Set<String> excludes =  new HashSet<>();
        Set<String> virtuals = new HashSet<>();
        excludes.add("className"); //always strip out class name
        DBObject output = new BasicDBObject();
        try {
            MappedClass mc = getMappedClass(Class.forName((String) data.get("className"), true, cl)); 
            for (MappedField field : mc.getMappedFields()){
                if (field.hasAnnotation(DBOnly.class) || field.hasAnnotation(Replace.class)){
                    excludes.add(field.getNameToStore());
                } else if (field.getField() == null && field.getNameToStore().startsWith(ClothoMappedField.VIRTUAL_FIELD_PREFIX)){
                    virtuals.add(field.getNameToStore());
                }
            }
        } catch (ClassNotFoundException ex) {
            logger.error("Schema not found", ex); 
        }
        
        for (Object key : data.keySet()){
            if (!excludes.contains(key)){
                String keyString = (String) key;
                if (virtuals.contains(keyString)){
                    keyString = keyString.substring(ClothoMappedField.VIRTUAL_FIELD_PREFIX.length());
                }
                output.put(keyString, data.get(key));
            }
        }

        return output;
    }
    
    public BSONObject toJSON(ObjBase obj){
        return toJSON(toDBObject(obj).toMap());
    }

    @Override
    public MappedClass getMappedClass(Object obj) {
        if (obj == null)
            return null;

        Class type = (obj instanceof Class) ? (Class) obj : obj.getClass();
        if (ProxyHelper.isProxy(obj))
            type = ProxyHelper.getReferentClass(obj);

        MappedClass mc = mappedClasses.get(type.getName());
        if (mc == null) {
            mc = new ClothoMappedClass(type, this);
            // no validation
            addMappedClass(mc, false);
        }
        return mc;
    }

    @Override
    public MappedClass addMappedClass(Class c) {
        MappedClass mc = new ClothoMappedClass(c, this);
        return addMappedClass(mc, true);
    }
}

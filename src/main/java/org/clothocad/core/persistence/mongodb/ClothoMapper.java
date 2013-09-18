/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.mongodb;

import com.github.jmkgreen.morphia.annotations.Reference;
import com.github.jmkgreen.morphia.mapping.DefaultMapper;
import com.github.jmkgreen.morphia.mapping.MappedClass;
import com.github.jmkgreen.morphia.mapping.MappedField;
import com.github.jmkgreen.morphia.mapping.MapperOptions;
import com.github.jmkgreen.morphia.mapping.cache.EntityCache;
import com.github.jmkgreen.morphia.mapping.lazy.proxy.ProxyHelper;
import com.mongodb.DB;
import com.mongodb.DBObject;
import com.mongodb.DBRef;
import java.util.List;
import java.util.Map;
import javax.inject.Inject;
import javax.persistence.EntityNotFoundException;
import lombok.Getter;
import lombok.Setter;
import org.bson.types.ObjectId;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.persistence.Add;
import org.clothocad.core.persistence.DBOnly;
import org.clothocad.core.persistence.Replace;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
/**
 *
 * @author spaige
 */
//TODO: preferentially remove schema classes from cache
import org.clothocad.core.aspects.JSONSerializer;
import org.clothocad.core.persistence.Adds;
import org.clothocad.core.persistence.DBClassLoader;
import org.clothocad.core.persistence.IdUtils;

public class ClothoMapper extends DefaultMapper implements JSONSerializer {

    protected BSONScrubber scrubber = new BSONScrubber(this);

    public ClothoMapper(MapperOptions opts) {
        super(opts);
    }
    
    @Inject
    public ClothoMapper(DBClassLoader cl){
       super(defaultOptions);
       this.cl = cl; 
    }
    
    static final Logger logger = LoggerFactory.getLogger(ClothoMapper.class);
    static final MapperOptions defaultOptions = new MapperOptions();

    static {
        defaultOptions.objectFactory = new PolymorphicObjectFactory();
    }

    static {
        //add our annotations to the interesting annotations lists
        MappedClass.interestingAnnotations.add(Add.class);
        MappedClass.interestingAnnotations.add(Adds.class);
        MappedField.interestingAnnotations.add(DBOnly.class);
        MappedField.interestingAnnotations.add(Replace.class);
    }
    @Getter
    @Setter
    //XXX: circular reference (this -> cl -> persistor -> this)
    protected ClassLoader cl;

    @Setter
    protected DB db;
    
    
    @Override
    public Map<String, Object> toJSON(Map data) {
        return scrubber.scrub(data);
    }

    @Override
    protected void readMappedField(DBObject dbObject, MappedField mf, Object entity, EntityCache cache) {
        if (mf.hasAnnotation(Replace.class)) {
            Replace annotation = mf.getAnnotation(Replace.class);
            try {
                entity.getClass().getMethod(annotation.decoder(), Map.class).invoke(entity, dbObject);
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
        } else {
            //massage any String references to DbRefs
            //TODO: Make this work with proxies
            if (mf.hasAnnotation(Reference.class)){
                Object fieldValue = dbObject.get(mf.getNameToStore());
                if (mf.isMultipleValues()){
                    //XXX: that second check shouldn't do anything
                    if (fieldValue != null){
                        List values = (List) fieldValue;
                        for (int i=0; i<values.size(); i++){
                            Object value = values.get(i);
                            if (value instanceof String){
                                values.set(i, convertToDBRef(value.toString()));
                            }
                        }
                    }
                }
                else {
                    if (!(fieldValue instanceof DBRef) && fieldValue != null){
                        if (fieldValue instanceof String){
                            dbObject.put(mf.getNameToStore(), convertToDBRef(fieldValue.toString()));
                        }
                    }
                }
            }
            
            super.readMappedField(dbObject, mf, entity, cache); //To change body of generated methods, choose Tools | Templates.
        }
    }

    protected DBRef convertToDBRef(String ref){
        try {
            return new DBRef(db, "data", new ObjectId(ref));
        } catch (IllegalArgumentException e){
            ObjectId id = IdUtils.resolveSelector(ref, false);
            if (id == null){
                throw new EntityNotFoundException("Could not find schema for selector "+ ref);
            }
            return new DBRef(db, "data", id);
        }
    }
    
    @Override
    public Map<String, Object> toJSON(ObjBase obj) {
        return toJSON(toDBObject(obj).toMap());
    }

    @Override
    public MappedClass getMappedClass(Object obj) {
        if (obj == null) {
            return null;
        }

        Class type = (obj instanceof Class) ? (Class) obj : obj.getClass();
        if (ProxyHelper.isProxy(obj)) {
            type = ProxyHelper.getReferentClass(obj);
        }

        MappedClass mc = mappedClasses.get(type.getName());
        if (mc == null) {
            logger.debug("Mapping class {}", type.getCanonicalName());
            mc = new ClothoMappedClass(type, this);
            // no validation
            addMappedClass(mc, false);
        }
        return mc;
    }

    @Override
    public MappedClass addMappedClass(Class c) {
        logger.debug("Mapping class {}", c.getCanonicalName());
        MappedClass mc = new ClothoMappedClass(c, this);
        return addMappedClass(mc, true);
    }

    @Override
    public List<Object> toJSON(List data) {
        return scrubber.scrub(data, null);
    }
}

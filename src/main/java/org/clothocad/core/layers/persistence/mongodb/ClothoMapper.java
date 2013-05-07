/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.layers.persistence.mongodb;

import com.github.jmkgreen.morphia.annotations.PrePersist;
import com.github.jmkgreen.morphia.annotations.PreSave;
import com.github.jmkgreen.morphia.mapping.DefaultMapper;
import com.github.jmkgreen.morphia.mapping.MappedClass;
import com.github.jmkgreen.morphia.mapping.MappedField;
import static com.github.jmkgreen.morphia.mapping.Mapper.CLASS_NAME_FIELDNAME;
import com.github.jmkgreen.morphia.mapping.MappingException;
import com.github.jmkgreen.morphia.mapping.cache.EntityCache;
import com.mongodb.BasicDBObject;
import com.mongodb.DBObject;
import java.util.Map;

/**
 *
 * @author spaige
 */
public class ClothoMapper extends DefaultMapper {

    @Override
    public Object fromDBObject(Class entityClass, DBObject dbObject, EntityCache cache) {
        return super.fromDBObject(entityClass, dbObject, cache); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public DBObject toDBObject(Object entity, Map<Object, DBObject> involvedObjects) {
        return toDBObject(entity, involvedObjects); //To change body of generated methods, choose Tools | Templates.
    }
    //TODO:
    //example of whole-object validation
    //transform and filter json
        //remove internal references (to class names)
        //normalize w/ sharable json
        //filter DB-specific info (like isDeleted, classData)
        
    
    DBObject toDBObject(Object entity, Map<Object, DBObject> involvedObjects, boolean lifecycle) {

        DBObject dbObject = new BasicDBObject();
        MappedClass mc = getMappedClass(entity);

        if (mc.getEntityAnnotation() == null || !mc.getEntityAnnotation().noClassnameStored())
            dbObject.put(CLASS_NAME_FIELDNAME, entity.getClass().getName());

        if (lifecycle)
            dbObject = (DBObject) mc.callLifecycleMethods(PrePersist.class, entity, dbObject, this);

        for (MappedField mf : mc.getMappedFields()) {
            try {
                writeMappedField(dbObject, mf, entity, involvedObjects);
            } catch (Exception e) {
                throw new MappingException("Error mapping field:" + mf.getFullName(), e);
            }
        }
        if (involvedObjects != null)
            involvedObjects.put(entity, dbObject);

        if (lifecycle)
            mc.callLifecycleMethods(PreSave.class, entity, dbObject, this);

        return dbObject;
    }
}

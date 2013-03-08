/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.layers.persistence.mongodb;

import com.github.jmkgreen.morphia.logging.Logr;
import com.github.jmkgreen.morphia.logging.MorphiaLoggerFactory;
import com.github.jmkgreen.morphia.mapping.DefaultCreator;
import com.github.jmkgreen.morphia.mapping.Mapper;
import com.mongodb.DBObject;

/**
 *
 * @author spaige
 */
public class PolymorphicObjectFactory extends DefaultCreator {
    private static final Logr log = MorphiaLoggerFactory.get(DefaultCreator.class);
    
    @Override 
    public Object createInstance(Class clazz, DBObject dbObj) {
        Class c = clazz;
        if (c == null)
            c = getClass(dbObj);
        return createInstance(c);
    }
    
    private Class getClass(DBObject dbObj) {
        // see if there is a className value
        String className = (String) dbObj.get(Mapper.CLASS_NAME_FIELDNAME);
        Class c = null;
        if (className != null) {
            // try to Class.forName(className) as defined in the dbObject first,
            // otherwise return the entityClass
            try {
                c = Class.forName(className, true, getClassLoaderForClass(className, dbObj));
            } catch (ClassNotFoundException e) {
                if (log.isWarningEnabled())
                    log.warning("Class not found defined in dbObj: ", e);
            }
        }
        return c;
    }
}

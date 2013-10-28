/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.mongodb;

import com.github.jmkgreen.morphia.logging.Logr;
import com.github.jmkgreen.morphia.logging.MorphiaLoggerFactory;
import com.github.jmkgreen.morphia.mapping.DefaultCreator;
import com.github.jmkgreen.morphia.mapping.Mapper;
import com.mongodb.DBObject;
import java.lang.reflect.Modifier;
import org.clothocad.core.persistence.IdUtils;
import org.clothocad.core.schema.BuiltInSchema;
import org.clothocad.core.schema.Schema;

/**
 *
 * @author spaige
 */
public class PolymorphicObjectFactory extends DefaultCreator {
    private static final Logr log = MorphiaLoggerFactory.get(PolymorphicObjectFactory.class);
    
    @Override 
    public Object createInstance(Class clazz, DBObject dbObj) {
        Class c = clazz;
        if (c == null  || Modifier.isAbstract(c.getModifiers())) //punt to stored class if clazz is abstract
            c = getClass(dbObj);
        return createInstance(c);
    }
    
    private Class getClass(DBObject dbObj) {

               
        // see if there is a className value
        String className = (String) dbObj.get(Mapper.CLASS_NAME_FIELDNAME);
        
        //is this describing a builtin?
        if (className.equals("BuiltInSchema")){
            /*try {
                return Class.forName(dbObj.get("c").toString());
            } catch (ClassNotFoundException ex) {
                log.warning("Built in class described by DB object not found: {}", dbObj.toString());
            }*/
            return BuiltInSchema.class;
        }
        
        Class c = null;
        if (className != null) {
            // try to Class.forName(className) as defined in the dbObject first,
            // otherwise return the entityClass
            try {
                c = Class.forName(className, true, getClassLoaderForClass(className, dbObj));
            } catch (ClassNotFoundException e) {
                //maybe it's just a schema name, try converting it
                try {
                    String convertedClassName = Schema.getBinaryName(IdUtils.resolveSelector(className, false));
                    //XXX: is still crashing, try getting a trail or something
                    c = Class.forName(convertedClassName, true, Schema.cl);
                } catch (ClassNotFoundException ex) {
                    if (log.isWarningEnabled())
                    log.warning("Class not found defined in dbObj: ", e);                   
                }
            }
        }
        return c;
    }
}

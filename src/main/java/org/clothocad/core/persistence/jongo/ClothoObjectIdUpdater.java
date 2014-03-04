/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jongo;

import org.bson.types.ObjectId;
import org.clothocad.core.datums.ObjBase;
import org.jongo.ReflectiveObjectIdUpdater;

/**
 *
 * @author spaige
 */
public class ClothoObjectIdUpdater extends ReflectiveObjectIdUpdater{

    public ClothoObjectIdUpdater(IdFieldSelector idFieldSelector) {
        super(idFieldSelector);
    }

    @Override
    public Object getId(Object pojo) {
        if (pojo instanceof ObjBase){
            ObjBase obj = (ObjBase) pojo;
            return obj.getId().toString();
        }
        return super.getId(pojo); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void setObjectId(Object newPojo, ObjectId id) {
        if (newPojo instanceof ObjBase){
            ObjBase obj = (ObjBase) newPojo;
            obj.setId(new org.clothocad.core.datums.ObjectId(id.toString()));
            return;
        }
        super.setObjectId(newPojo, id); //To change body of generated methods, choose Tools | Templates.
    }
    
    
}

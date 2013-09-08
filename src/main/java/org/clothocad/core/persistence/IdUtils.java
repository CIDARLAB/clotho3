/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence;

import javax.inject.Inject;
import org.bson.types.ObjectId;
import org.clothocad.core.datums.ObjBase;

/**
 * lame hack to handle the need to look up various ids in the data layer
 * persistor injected by clotho module
 * @author spaige
 */
public class IdUtils {
    @Inject
    static Persistor persistor;
    
    public static ObjBase get(ObjectId id){
        return persistor.get(ObjBase.class, id);
    }
}

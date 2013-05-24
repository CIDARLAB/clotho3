/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.aspects;

import java.util.Map;
import org.bson.BSONObject;
import org.clothocad.core.datums.ObjBase;

/**
 *
 * @author spaige
 */
public interface JSONSerializer {
    public BSONObject toJSON(ObjBase obj);
    public BSONObject toJSON(Map data);
}

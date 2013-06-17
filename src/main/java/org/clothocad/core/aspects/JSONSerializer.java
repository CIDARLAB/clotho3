/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.aspects;

import java.util.Map;
import org.clothocad.core.datums.ObjBase;

/**
 *
 * @author spaige
 */
public interface JSONSerializer {
    public Map<String,Object> toJSON(ObjBase obj);
    public Map<String,Object> toJSON(Map data);
}

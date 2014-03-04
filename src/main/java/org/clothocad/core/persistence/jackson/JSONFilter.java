/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.jackson;

import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonUnwrapped;
import java.util.HashMap;
import java.util.Map;
import org.clothocad.core.datums.ObjBase;

/**
 *
 * @author spaige
 */
public class JSONFilter extends ObjBase{
    
    
    
    @JsonUnwrapped
    private Map<String,Object> _fields = new HashMap<>();
    
    @JsonAnySetter
    public void stashField(String name, Object o){
        if (name.equals("_id")) name = "id";
        
        
        _fields.put(name, o);
    }
    
}

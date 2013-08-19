/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.model;

import java.util.List;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.AllArgsConstructor;
import org.clothocad.core.aspects.Proctor.Module;
import org.clothocad.core.datums.SharableObjBase;

/**
 *
 * @author jcanderson
 */
@NoArgsConstructor
@AllArgsConstructor
public class Trail extends SharableObjBase {
    
    public Trail(String name, String description, List<Module> contents){
        super(name, null, description);
        this.contents = contents;
    }
    
    @Getter
    private List<Module> contents;
}

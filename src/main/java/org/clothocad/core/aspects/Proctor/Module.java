/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.aspects.Proctor;

/**
 *
 * @author jcanderson
 */

import java.util.List;
import lombok.AllArgsConstructor;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

@AllArgsConstructor
@NoArgsConstructor
public class Module {
    @Setter
    @Getter
    private String module_title;
    @Getter
    private List<Paver> pavers;
    
}

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.model;

import java.util.List;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.AllArgsConstructor;
import lombok.Setter;
import org.clothocad.core.aspects.Proctor.Module;
import org.clothocad.core.datums.ObjBase;
import org.json.JSONObject;

/**
 *org.clothocad.model.Trail
 * @author jcanderson
 */
@NoArgsConstructor
@AllArgsConstructor
public class Trail extends ObjBase {

	/**
    public Trail(String nameOrTitle, String description, List<Module> contents) {
        super(nameOrTitle, null);
        this.description = description;
        this.title = nameOrTitle;
        this.contents = contents;
    }
    **/

    @Setter
    @Getter
    private String title;
    @Setter
    @Getter
    private String description;
    @Getter
    private List<Module> contents;
}

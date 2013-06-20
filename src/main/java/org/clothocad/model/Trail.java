/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.model;

import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.AllArgsConstructor;
import lombok.Setter;
import org.clothocad.core.aspects.Proctor.Module;
import org.clothocad.core.datums.Sharable;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;

/**
 *org.clothocad.model.Trail
 * @author jcanderson
 */
@NoArgsConstructor
@AllArgsConstructor
public class Trail extends Sharable {

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

    @Override
    public boolean validate(JSONObject obj) {
        return true;    
    }
}

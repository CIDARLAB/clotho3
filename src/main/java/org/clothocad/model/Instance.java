/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.model;

import lombok.Getter;
import lombok.NoArgsConstructor;
import org.clothocad.core.datums.Sharable;
import org.json.JSONObject;

/**
 *org.clothocad.model.Instance
 * @author jcanderson
 */
@NoArgsConstructor
public class Instance extends Sharable {

    public Instance(JSONObject data) {
        json = data;
    }
    @Getter
    private JSONObject json;
    
    @Override
    public boolean validate(JSONObject obj) {
        //run a scriptengine on the data and find out
        return true;
    }
}

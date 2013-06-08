/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.model;


import java.util.logging.Level;
import java.util.logging.Logger;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.clothocad.core.datums.Sharable;
import org.json.JSONException;
import org.json.JSONObject;

/**
 *
 * @author jcanderson
 */
@NoArgsConstructor
public class Trail extends Sharable {

    /**
     * Constructor from raw data
     * @param data
     */
    
    @Setter
    @Getter
    private String data;
    
    //TODO:unique name constraint
    public Trail( String name, JSONObject data ) {
        super(name, null);
        this.data = data.toString();
    }
    
    @Override
    public JSONObject toJSON() {
        try {
            return new JSONObject(data);
        } catch (JSONException ex) {
            ex.printStackTrace();
            return null;
        }
    }

}

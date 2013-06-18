/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.model;

import java.util.logging.Level;
import java.util.logging.Logger;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.AllArgsConstructor;
import lombok.Setter;
import org.clothocad.core.datums.Sharable;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;

/**
 *
 * @author jcanderson
 */
@NoArgsConstructor
public class Trail extends Sharable {
    
    public Trail(String nameOrTitle, String description, JSONArray contents) {
        super(nameOrTitle, null);
        this.description = description;
        this.title = nameOrTitle;
        this.contents = contents.toString();
    }

    @Setter
    @Getter
    private String title;
    @Setter
    @Getter
    private String description;

    private String contents;
    
    public JSONArray getContents() {
        try {
            return new JSONArray(contents);
        } catch (JSONException ex) {
            return null;
        }
    }
    
    public void setContents(JSONArray data) {
        contents = data.toString();
    }
    
    @Override
    public boolean validate(JSONObject obj) {
        return true;    
    }
}

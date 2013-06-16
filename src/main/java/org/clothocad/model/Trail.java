/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.model;

import java.util.ArrayList;
import java.util.List;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.AllArgsConstructor;
import lombok.Setter;
import org.clothocad.core.aspects.Proctor.Paver;
import org.clothocad.core.datums.Sharable;
import org.json.JSONException;
import org.json.JSONObject;

/**
 *
 * @author jcanderson
 */
@NoArgsConstructor
@AllArgsConstructor
public class Trail extends Sharable {

    /**
     * Constructor from raw data
     * @param data
     */
    /** JSON of a Trail looks like:
    "uuid" : "trail_123",
    "author" : "UC Berkeley",       
    
    "title" : "Biosafety Module",
    "description" : "desc",
    "contents" : []
    }
     */
    

    @Setter
    @Getter
    private String title;
    @Setter
    @Getter
    private String description;
    @Setter
    @Getter
    private List<Content> contents = new ArrayList<Content>();
    
    @AllArgsConstructor
    public static class Content {
        @Setter
        @Getter
        private String module_title;
        @Setter
        @Getter
        private List<Paver> pavers = new ArrayList<Paver>();
    }


    @Override
    public boolean validate(JSONObject obj) {
        return true;    
    }
    

}

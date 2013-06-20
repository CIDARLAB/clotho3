/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.aspects.Proctor;

import java.util.logging.Level;
import java.util.logging.Logger;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.json.JSONException;
import org.json.JSONObject;

/**
 *
 * @author jcanderson
 */
@NoArgsConstructor
public class TemplatePaver extends Paver {

    public TemplatePaver(String paver_title, String template) {
        super(paver_title, "template");
        this.template = template;
        testObj = new JSONObject();
        try {
            testObj.put("testy", template);
        } catch (JSONException ex) {
            Logger.getLogger(TemplatePaver.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    @Getter
    @Setter
    private String template;
    
    @Getter
    @Setter
    private JSONObject testObj;
}

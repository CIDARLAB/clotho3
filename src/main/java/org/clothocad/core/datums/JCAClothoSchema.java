/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.datums;

import org.json.JSONObject;

/**
 *
 * @author jcanderson
 */
public interface JCAClothoSchema {
    public boolean validate(JSONObject obj);
}

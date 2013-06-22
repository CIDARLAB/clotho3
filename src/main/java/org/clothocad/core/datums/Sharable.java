/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.datums;

import org.clothocad.model.Person;
import org.json.JSONObject;

/**
 *
 * @author jcanderson
 */
public interface Sharable  {

    //Functional requirements of all Sharables
    public JSONObject toJSON();
    public String getId();
    
    //Metadata for all Sharables
    Person getAuthor();
    public String getIcon();
    public String getName();
    public String getDescfription();
}

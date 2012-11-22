/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.datums;

import org.clothocad.core.datums.objbases.Person;
import org.json.JSONObject;

/**
 *
 * @author jcanderson
 */
public interface Sharable extends Datum {
    JSONObject toJSON();
    Person extractAuthor();
    //NEED BACK THE SHARING RULES FROM OLD DATUM
    
    public SharableType type();

    /**
     * Each Sharable class is responsible for handling its own permissions
     * 
     * @param newvalue
     * @param requestor
     * @return 
     */
    public boolean set(JSONObject newvalue, Person requestor, Doo doo);

    public static enum SharableType {SCHEMA, VIEW, FUNCTION, INSTANCE, BADGE, QUIZ, TRAIL, CONTENT};
}

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.datums;

import java.util.Iterator;
import lombok.Getter;
import org.clothocad.core.aspects.Persistor;
import org.clothocad.model.Person;
import org.json.JSONObject;

/**
 *
 * @author jcanderson
 */
public abstract class Sharable extends ObjBase {

    @Getter
    private Person author;
        
    public Sharable() {}
	
    public Sharable(String name, Person author) {
        super(name);
        this.author = author;
    }
    
    @Override
    public String toString() {
        JSONObject obj = this.toJSON();
        return obj.toString();
    }
    
/**
     * Stephanie, I need you to implement this in various spots
     * @param obj
     * @return 
     */
    public abstract boolean validate(JSONObject obj);
}

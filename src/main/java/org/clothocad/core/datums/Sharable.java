/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.datums;

import lombok.Getter;
import org.clothocad.model.Person;

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
    
}

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence;

import lombok.Getter;
import lombok.Setter;

import org.clothocad.model.Format;
import org.clothocad.model.Part;
import org.clothocad.model.Person;
import org.clothocad.model.Sequence;

/**
 *
 * @author spaige
 */
public class ExtendedPart extends Part {
    //basic constructor 
    //Part generateBasic(String name, String shortdescription, String seq, Format form, Person author) 
    //    private Part(String name, String shortDescription, String seq, Format form, Person author) {
    public ExtendedPart(String name, String desc, Sequence seq, Format form, Person author){
        super(name, desc, seq, author);
        setFormat(form);
    }
    
    //composite constructor
    
    
    @Getter
    @Setter
    private String additionalParameters;
    
}

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.model;

import java.util.List;

import lombok.NoArgsConstructor;

import org.clothocad.model.Format;
import org.clothocad.model.Part;

/**
 *
 * @author spaige
 * @author Nicholas Roehner
 */
@NoArgsConstructor
public class FreeForm implements Format {

	@Override
    public boolean checkPart(Part p) {
        return true;
    }

	@Override
    public boolean checkComposite(List<Part> composition) {
        return true;
    }

	@Override
    public Part generateCompositePart(String name, List<Part> composition, Person author) {
        StringBuilder builder = new StringBuilder();
        for (Part part : composition){
            builder.append(part.getSequence().getSequence());
        }
        Part compositePart = new Part(name, new SimpleSequence(builder.toString(), author), author);
        compositePart.setFormat(this);
        Assembly assembly = new Assembly(compositePart);
        for (Part subPart : composition) {
        	assembly.getParts().add(subPart);
        }
        return compositePart;
    }
	
	@Override
    public Part generateCompositePart(String name, String description, List<Part> composition, Person author) {
        StringBuilder builder = new StringBuilder();
        for (Part part : composition){
            builder.append(part.getSequence().getSequence());
        }
        Part compositePart = new Part(name, description, new SimpleSequence(builder.toString(), author), author);
        compositePart.setFormat(this);
        Assembly assembly = new Assembly(compositePart);
        for (Part subPart : composition) {
        	assembly.getParts().add(subPart);
        }
        return compositePart;
    }
    
}

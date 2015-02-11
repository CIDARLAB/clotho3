/*
Copyright (c) 2009 The Regents of the University of California.
All rights reserved.
Permission is hereby granted, without written agreement and without
license or royalty fees, to use, copy, modify, and distribute this
software and its documentation for any purpose, provided that the above
copyright notice and the following two paragraphs appear in all copies
of this software.

IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.
createQuery
THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
ENHANCEMENTS, OR MODIFICATIONS..
 */
package org.clothocad.model;

import java.util.LinkedList;
import java.util.List;

import javax.validation.constraints.NotNull;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.model.Feature.FeatureRole;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

/**
 * @author jcanderson
 * @author Nicholas Roehner
 */
@NoArgsConstructor
public class Part extends SharableObjBase {
	
	@Getter
	@NotNull
	private Format format;
	
	@Getter
	@Setter
	private List<Assembly> assemblies = new LinkedList<Assembly>();
	
	@NotNull
	@Getter
	@Setter
	private Sequence sequence;
	
	@Getter
    @Setter
    private boolean isForwardOrientation;
	
	@Getter
	@Setter
	private Part parentPart;
    
    public Part(String name, String description, Sequence sequence, Person author){
        super(name, author, description); 
        this.sequence = sequence;
    }
    
    public Part(String name, Sequence sequence, Person author){
        super(name, author); 
        this.sequence = sequence;
    }

    /**
     * Change the Format of the Part
     * @param format new Format for the Part
     */
    public void setFormat(Format format) {
        if (format != null && format.checkPart(this)) {
        	this.format = format;
        }
    }
    
    public List<FeatureRole> getRoles() {
    	List<Annotation> annotations = sequence.getAnnotations();
    	List<FeatureRole> roles = new LinkedList<FeatureRole>();
    	for (Annotation annotation : annotations) {
    		Feature feature = annotation.getFeature();
    		if (feature != null) {
    			roles.add(feature.getRole());
    		}
    	}
    	return roles;
    }
    
}

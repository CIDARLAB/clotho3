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

import java.util.List;

import javax.validation.constraints.AssertTrue;
import javax.validation.constraints.NotNull;

import org.clothocad.core.datums.SharableObjBase;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import lombok.extern.slf4j.Slf4j;

/**
 * @author jcanderson
 */
@NoArgsConstructor
@Slf4j
public abstract class Part extends SharableObjBase {
	
	@Getter
	@NotNull
	private Format format;

	@Getter
	@Setter
	@Deprecated
	// Type (actually role) should be derived from Features annotating the Part's Sequence.
	// Part should only store data on composition of sequences.
	private PartFunction type;

	@Getter
	@Setter
	@Deprecated
	// Risk should be derived from the Part's Sequence.
	// Sequences and their Annotations are connection between
	// genetic structure and function.
	private short riskGroup;
    
    protected Part(String name, String description, Format format, Person author){
        super(name, author, description); 
        this.format = format;
    }
    
    /**
     * Call this method to construct a new basic Part. It will check that the Part obeys its Format. 
     * If not, then the method returns null.
     *
     * @param name nickname of Part, for example "R0040"
     * @param description description of Part, for example "TetR-repressible promoter"
     * @param sequence sequence of the Part, for example "cgaaggcaggacacacatg"
     * @param format Format of the Part
     * @param author author of the Part, for example Person with name "Bob"
     */
    public static Part generateBasic(String name, String description, String sequence, Format format, Person author) {
        return new BasicPart(name, description, sequence, format, author);
    }

    /**
     * Create composite Part from Part list
     *
     * @param composition list of Parts to compose
     * @param name nickname of Part, for example "C0040"
     * @param description description of Part, for example "coding sequence for TetR repressor with LVA tail"
     * @param sequence sequence of the Part, for example "cgaaggcaggacacacatg"
     * @param format Format of the Part
     * @param author author of the Part, for example Person with name "Bob"
     */
    public static Part generateComposite(List<Part> composition, Format format, Person author, String name, String description) {
        return new CompositePart(composition, format, author, name, description);
    }
   
    @AssertTrue
    public abstract boolean checkFormat();
    
    public String getFormatName(){
        return format.getClass().getSimpleName();
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

    public abstract NucSeq getSequence();

    public static Part retrieveByName(String name) {
        //query connection for one part whose name contains the provided string    
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @AssertTrue
    private boolean checkDBConstraints(){
        //name must be unique
        //format(by name) + sequence(by content)combination should be unique
        return true;
    }
    
    @Deprecated
    // Type should be derived from Features annotating the Part's Sequence.
 	// Part should only store data on composition of sequences.
    public static enum PartFunction {
        //XXX: composite is not a function
        CDS, RBS, PROMOTER, TERMINATOR, COMPOSITE;
    }
}

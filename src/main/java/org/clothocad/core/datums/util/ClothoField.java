/*
 * 
Copyright (c) 2010 The Regents of the University of California.
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

THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
ENHANCEMENTS, OR MODIFICATIONS..
 */

package org.clothocad.core.datums.util;

import lombok.Getter;
import org.clothocad.core.datums.Function;
/**
 * @author John Christopher Anderson
 */

@Getter
public class ClothoField {
    
    private ClothoField() {}
    
    public ClothoField(String name, Class type, String example, String description, Function validate, boolean reference, int access) {
        this.name = name;
        this.type = type;
        this.example = example;
        this.access = access; 
        this.reference = reference;
        this.validate = validate;
        this.description = description;
    }
    

    private String name;
    private Class type;
    private String example;   //A string representation/explanation of an expected value
    private int access;  //uses asm opcodes
    private boolean reference;
    private Function validate;
    
    //metadata
    private String description;
}
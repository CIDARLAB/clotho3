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

package org.clothocad.core.testers;

import javax.script.ScriptException;
import org.clothocad.core.datums.Doo;
import org.clothocad.core.layers.communication.mind.Mind;

/**
 * @author John Christopher Anderson
 */


public class MindTest {
    
    public static void main(String[] args) throws ScriptException {
        Mind mind = new Mind();
        /* TODO: put Doo back */
        /* Doo doo = new Doo(null, false); */
        
        mind.runCommand("clotho.say('hi');");
        mind.runCommand("clotho.say('makeabuddy');");
        mind.runCommand("var asst = clotho.get('makeabuddyassistantisuuid');");
        mind.runCommand("clotho.say('I have assistant: ' + asst);");
        
        mind.runCommand("var cind = clotho.get('cindysuobjbasetypeperson');");
        mind.runCommand("clotho.say('I have cind: ' + cind);");
        
        mind.runCommand("clotho.run(asst, 'cindysuobjbasetypeperson');");

        mind.runCommand("var dd = clotho.create('personschemaisuuid', '{\"name\":\"dougd\",\"email\":\"doug@bob.com\"}');");
        mind.runCommand("clotho.say('I have dd: ' + dd);");
        mind.runCommand("clotho.run(asst, dd);");
    }
}

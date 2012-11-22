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

package Testers;

import javax.script.ScriptException;
import org.clothocad.core.aspects.Communicator.Mind.Mind;
import org.clothocad.core.datums.Doo;

/**
 * @author John Christopher Anderson
 */


public class MindTest {
    
    public static void main(String[] args) throws ScriptException {
        Mind mind = new Mind();
        /* TODO: put Doo back */
        /* Doo doo = new Doo(null, false); */
        
        mind.runCommand(null, "clotho.say('hi');");
        mind.runCommand(null, "clotho.say('makeabuddy');");
        mind.runCommand(null, "var asst = clotho.get('makeabuddyassistantisuuid');");
        mind.runCommand(null, "clotho.say('I have assistant: ' + asst);");
        
        mind.runCommand(null, "var cind = clotho.get('cindysuobjbasetypeperson');");
        mind.runCommand(null, "clotho.say('I have cind: ' + cind);");
        
        mind.runCommand(null, "clotho.run(asst, 'cindysuobjbasetypeperson');");

        mind.runCommand(null, "var dd = clotho.create('personschemaisuuid', '{\"name\":\"dougd\",\"email\":\"doug@bob.com\"}');");
        mind.runCommand(null, "clotho.say('I have dd: ' + dd);");
        mind.runCommand(null, "clotho.run(asst, dd);");
    }
}

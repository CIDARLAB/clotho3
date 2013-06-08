/*
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
ENHANCEMENTS, OR MODIFICATIONS.
 */

package org.clothocad.core.testers;

import org.clothocad.core.aspects.Persistor;
import org.clothocad.core.datums.Datum;
import org.clothocad.core.layers.communication.mind.Mind;

/**
 * Tests the NameSpace functionality of the Mind
 * Create a Mind, populate it with stuff, clear it, then run something namespaced
 */
public class T6 {
    /*private static final int counter = 2;
    
    public static void main(String[] args) {
        //Create or get a static mind reference
        Mind mind = null;
        Datum dat = Collector.get().getDatum("mind-instance-test7-uuid");
        if(dat!=null) {
            mind = (Mind) dat;
        } else {
            mind = new Mind();
            mind.SUPERILLEGAL_SETUUID("mind-instance-test7-uuid");
            mind.save();
        }

        /* TODO: put Doo back */
        /* Doo doo = new Doo(null, false); */
        /*
        //Put in some commands
        Datum personSchema = Person.getSchema() ;
        for(int i=0; i<counter; i++) {
            mind.runCommand(null, "var person" + i + " = clotho.get('" + personSchema.getId() + "');");
        }
        
        Datum Cindy = Collector.get().getDatum("specific-cindysu-is-uuid");
        for(int i=0; i<counter; i++) {
            mind.runCommand(null, "var cindy" + i + " = clotho.get('" + Cindy.getId() + "');");
        }
        
        //Let's see what the engine has for a recently-created object
        System.out.println("T7 testing a command while js is still intact");
        mind.runCommand(null, "clotho.log(\"INFO\", 'T7: person1 is: ');");
        mind.runCommand(null, "clotho.log(\"INFO\", person1);");
        mind.runCommand(null, "var acopy = person1;");
        mind.runCommand(null, "clotho.log(\"INFO\", acopy);");
        
        //Clear the engine, which won't clear the namespace
        mind.clear();
        
        //See what now executes (just the namespaced tokens)
        Logger.log(Logger.Level.INFO, "T7 testing a command after cleared");
        mind.runCommand(null, "clotho.log(\"INFO\", 'cindy is still in namespace: ' + cindy1);");
        mind.runCommand(null, "clotho.log(\"INFO\", 'the copy generates an error now ' + acopy);");
        
        //Let's try it on a function:
        mind.runCommand(null, "var afunc = function(person) { clotho.log(\"INFO\", 'afunc says ' + person); };");
        mind.runCommand(null, "afunc(person1);");
        mind.clear();
        mind.runCommand(null, "afunc(person1);");
    }*/
}

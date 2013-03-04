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
ENHANCEMENTS, OR MODIFICATIONS..
 */

package org.clothocad.core.testers;

import org.clothocad.core.aspects.Collector;
import org.clothocad.core.datums.Instance;
import org.clothocad.core.datums.Schema;
import org.clothocad.core.util.Logger;

/**
 * Confirms that the Persistor finds the Feature and gfp datums
 * @author John Christopher Anderson
 */
public class T2 {
    public static void main(String[] args) {
        //Get the "Person" schema's contents
        Schema feature = (Schema) Collector.get().getDatum("specific-simplefeature-is-uuid");
        if (feature == null) {
            Logger.log(Logger.Level.FATAL, "Could not retrieve 'specific-simplefeature-is-uuid'");
            return;
        }
        Logger.log(Logger.Level.INFO, feature.getId());
        
        //Retrieve an instance of Person
        Instance gfp = (Instance) Collector.get().getDatum("specific-gfpuv-is-uuid");
        if (gfp == null) {
            Logger.log(Logger.Level.FATAL, "Could not retrieve 'specific-gfpuv-is-uuid'");
            return;
        }
        Logger.log(Logger.Level.INFO, gfp.getId());
        System.out.println("T2:  All is good!");
    }
}

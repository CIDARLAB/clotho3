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

package org.clothocad.core.aspects;

/**
 * The Authenticator is responsible for keeping track of the Users and Persons
 * that are registered on this instance of Clotho.
 * 
 * It is responsible for mapping a person trying to log in to a user with a login
 * and password.
 * 
 * @author John Christopher Anderson
 */


public class Authenticator implements Aspect {
    
    public boolean isOnDomain(String personUUID) {
        return true;
    }

    public boolean authenticate(String personUUID, String shaipass) {
        //LOOKUP THE USER ASSOCIATED WITH THAT PERSONUUID BY HASHMAP
        //IF THE PASSWORD MATCHES THAT OF THE USER, RETURN TRUE, OTHERWISE FALSE
        return true;
    }
   
    /**
     * SINGLETON BITS
     * modeled after http://www.javaworld.com/javaworld/jw-04-2003/jw-0425-designpatterns.html
     */
    private Authenticator() { }
    
    public static Authenticator get() {
        return singleton;
    }

    private static final Authenticator singleton = new Authenticator();
}

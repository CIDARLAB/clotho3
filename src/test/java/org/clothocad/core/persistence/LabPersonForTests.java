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

THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
ENHANCEMENTS, OR MODIFICATIONS..
 */
package org.clothocad.core.persistence;

import org.clothocad.core.persistence.annotations.Reference;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.model.*;

/**
 *
 * @author J. Christopher Anderson
 */
@NoArgsConstructor
public class LabPersonForTests extends Person {

    @Getter
    @Setter
    @Reference
    private Lab lab;
    @Getter
    @Setter
    @Reference
    private Collection myCollection;
    
    @Getter
    @Setter
    private String givenName, surName, nickName, registryName, emailAddress, snailMailAddress;
    

    /**Constructor from raw data
     *
     * @param displayname = String of author such as "JCAnderson"
     * @param affiliation = String of affiliation such as "UC Berkeley"
     */
    //unique name criterion
    //valid or nonexistent email
    public LabPersonForTests( String displayname, Lab alab, String rawPassword ) {
        //XXX:  Do people have authors?
        super(displayname, null);
        lab = alab;
        //changePassword( rawPassword );
        myCollection = new Collection();
        //biography = new WikiText("");
    }

    





    /* GETTERS
     * */

    /**
     * Retrieve a Person ObjBase from the database using their display name.
     * This will only query Persons saved to the database, not local ones.
     * @param name
     * @return
     */
    public static Person retrieveByName( String name ) {
        throw new UnsupportedOperationException();
    }

    /**
     * It checks a password to see if it matches the user's password and returns
     * true if they match, otherwise returns false.
     * It only stores the encrypted version of the password using SHA-1 hashing.
     *
     * @param raw the raw password supplied by user
     * @return true if it's a match
     */


    /*private boolean hasChangeClearance() {
        //If it's a new Person, changes are OK
        if(_isBrandNew) {
            return true;
        }
        Person user = Collector.getCurrentUser();
        //If nobody's logged in, they aren't allowed to change things
        if(user==null) {
            JOptionPane.showMessageDialog( null, "You aren't logged in.  You aren't allowed to change this data.", "Forbidden data change", JOptionPane.OK_OPTION );
            return false;
        }
        //If a person is editing their own fields, that's ok
        if(user.getUUID().equals(this.getUUID())) {
            return true;
        }
        //If the logged in user is an admin, they can edit fields
        if(user._personDatum._isAdministrator) {
            return true;
        }
        JOptionPane.showMessageDialog( null, "Only the user herself or admins can change this data.", "Forbidden data change", JOptionPane.OK_OPTION );
        return false;
    }*/

    /**
     * Is the person an administrator?
     * @return true if they are
     */
   // public boolean isAdmin() {
   //     return _personDatum._isAdministrator;
   // }

    /**
     * Get the person's login name
     * @return a String
     */
    public String getDisplayName() {
        return getName();
    }







    /**
     * Get the personal Collection of this object
     * @return a Collection ObjBase
     */
    public Collection getHerCollection() {
        return myCollection;
    }


    @Override
    public String toString() {
        return getName();
    }
}

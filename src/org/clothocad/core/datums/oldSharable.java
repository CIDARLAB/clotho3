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

package org.clothocad.core.datums;

import org.json.JSONObject;

/**
 *  I considered pulling stuff from second life, but this really isn't that scenario
 * So, I just thought through how it has to be.  It's probably not right.
 * 
 * So, I'm not sure that all this sharing stuff is the same for all the different
 * Datums -- like the logic here might get a little complicated, so the interface
 * just requires that things provide getters for what the permissions status is
 * and a single setter that can change those values.  I don't think this interface
 * will ever have to change since it's just abstracting the security process not
 * encoding it.  So, hopefully using this interface will give us the flexibility to
 * just run with whatever we need, and the guis and such built atop this won't be 
 * broken if we entirely change how the security is best implemented.
 * 
 * @author jcanderson
 */

public interface oldSharable {
    /* Visible/Invisible to other users. (e.g. in search results)
     * For an Assistant, there will also be a canRun() method.
     */
    boolean canView(String userId);

    /* Dynamic content can/cannot be edited by person.
     * Changing UUIDs and such is never allowed from within the container.
     */
    boolean canEdit(String userId);

    /* Can/Cannot be duplicated on another domain?
     * (can include canEdit and canView rules within it)
     */
    boolean canShare(String userId);

    /* Can/Cannot sell object */
    boolean canSell(String userId);
    
    Instance getOwner();
    
    JSONObject getJSON();

    void set(JSONObject permissions, String userId);
    
    SharableType getSharableType();
    
    public enum SharableType {BADGE, VIEW, ASSISTANT, SCHEMA, OBJBASE, TRAIL}
}

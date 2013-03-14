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

package org.clothocad.core.aspects.Ambassador;

import java.util.UUID;
import org.clothocad.core.datums.Datum;
import org.clothocad.core.datums.ObjBase;

/**
 * A Sister is the representation of another ClothoCore
 * @author John Christopher Anderson
 */


public class Sister 
		extends ObjBase {
    public Sister(String id, String url, int port) {
        this.id = id;
        this.url = url;
        this.port = port;
    }

    public String getId() {
        return id;
    }

    public int getPort() {
        return port;
    }

    public TrustLevel getTrust() {
        return trust;
    }

    public String getUrl() {
        return url;
    }

    public String getSecurityCode() {
        return securityCode;
    }

    public void setPort(int port) {
        this.port = port;
    }

    public void setTrust(TrustLevel trust) {
        this.trust = trust;
    }

    public void setUrl(String url) {
        this.url = url;
    }
    
    public static enum TrustLevel {
        BLOCKED,  //Reject all requests from this ClothoCore
        UNKNOWN,  //I'm seeing this ClothoCore for the first time
        TRUSTED   //I've gone through the process of granting this Clotho full access
    }
    
    private final String id;  //the uuid of that ClothoCore in all of Clotho land
    private String url; //the current url of that ClothoCore
    private int port;   //the port that this ClothoCore uses for communication
    private TrustLevel trust = TrustLevel.UNKNOWN;
    private String securityCode = UUID.randomUUID().toString();
}

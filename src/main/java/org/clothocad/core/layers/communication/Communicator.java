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
package org.clothocad.core.layers.communication;

import java.util.Map;
import java.util.HashMap;
import org.clothocad.core.aspects.Aspect;
import org.clothocad.core.layers.communication.mind.Mind;
import org.json.JSONException;
import org.json.JSONObject;


/**
/**
 * The Connector and Communicator have somewhat similar functions but not identical.
 * 
 * Communicator is for speaking to a client device, like an ipad, smartphone, or browser.
 * As such, it involves a single user at a time and the communication is all the API stuff
 * but also a bunch of visual stuff like GUI tracking.  It is a high-resolution channel
 * of communication which assumes that both the client and server have a fully-entrusted
 * relationship.
 * 
 * Ambassador is for speaking server to server or clothocore to clothocore.  As such, it is
 * just for CRUD operations, relaying of search requests, and confirmation of badges.
 * So, it doesn't involve GUIs.  The critical difference is that the Connector does not assume
 * its partner is trustworthy.
 * 
 * 
 * @author John Christopher Anderson
 */

public final class Communicator 
		implements Aspect {
	
    public void sendClientMessage(String connection_id,
                                  String channel,
                                  String message) {
    	JSONObject json = new JSONObject();
    	try {
			json.put("message", message);
			
			// TODO:
	        //Router.get().sendMessage(connection_id, channel, json);
		} catch (JSONException e) {
			e.printStackTrace();
		}
    }
    
    public Mind getMind(String auth_key) {
        /* assumes auth_key has been validated */
        Mind mind = null;
        if (auth_minds.containsKey(auth_key)) {
            mind = auth_minds.get(auth_key);
        } else {
            /* make a new mind */
            mind = new Mind();
            auth_minds.put(auth_key, mind);
        }
        return mind;
    }


	/** DOUBLE-CHECKED LOCKING **/
	private static volatile Communicator communicator;

    public static Communicator get() {
    	Communicator c = communicator;
		if(c == null) {
			synchronized(Communicator.class) {
				c = communicator;
				if(c == null) {
					communicator = c = new Communicator();
				}
			}
		}
		return c;
    }

    /* maps from "auth_key" to Mind */
    private Map<String, Mind> auth_minds = new HashMap<String, Mind>();
    private UpdateIndex updateReg = new UpdateIndex();
}

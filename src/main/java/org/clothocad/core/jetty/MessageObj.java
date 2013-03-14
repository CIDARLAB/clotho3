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
package org.clothocad.core.jetty;

import org.json.JSONObject;
import org.json.JSONException;

import org.clothocad.core.util.Logger;

/**
 * Helper object to parse/serialize (socket_id, channel, message, flags)
 * 4-tuple into/from JSON.
 *
 * @author Kelvin Li
 */
class MessageObj {
    static MessageObj fromJSONString(String json_string) {
        MessageObj mobj = new MessageObj();
        try {
            JSONObject json = new JSONObject(json_string);
            mobj.socket_id = (String) json.get("socket_id");
            mobj.channel = (String) json.get("channel");
            mobj.message = (String) json.get("message");
            mobj.flags = (String) json.get("flags");
        } catch (JSONException e) {
            Logger.log(Logger.Level.WARN,
                       "failed parsing " + json_string,
                       e);
            return null;
        } catch (ClassCastException e) {
            Logger.log(Logger.Level.WARN,
                       "failed parsing " + json_string,
                       e);
            return null;
        }
        return mobj;
    }

    static MessageObj fromComponents(String socket_id,
                                     String channel,
                                     String message,
                                     String flags) {
        MessageObj mobj = new MessageObj();
        mobj.socket_id = socket_id;
        mobj.channel = channel;
        mobj.message = message;
        mobj.flags = flags;
        return mobj;
    }

    JSONObject getJSON() {
        try {
            JSONObject json = new JSONObject();
            json.put("socket_id", socket_id);
            json.put("channel", channel);
            json.put("message", message);
            json.put("flags", flags);
            return json;
        } catch (JSONException e) {
            Logger.log(Logger.Level.WARN,
                       "failed constructing json",
                       e);
            return null;
        }
    }

    String getSocketID() {return socket_id;}
    String getChannel() {return channel;}
    String getMessage() {return message;}
    String getFlags() {return flags;}

    private String socket_id;
    private String channel;
    private String message;
    private String flags;
}

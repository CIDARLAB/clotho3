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

import java.util.Map;
import java.util.HashMap;
import java.io.IOException;
import org.eclipse.jetty.websocket.WebSocket;
import org.clothocad.core.util.Logger;

/**
 * @author Kelvin Li
 */
public class WSObj implements WebSocket.OnTextMessage {
    private WSObj() {}

    @Override
    public void onMessage(String data) {
    	/***
        MessageObj mobj = MessageObj.fromJSONString(data);
        if (Router.get().getInternalChannel().equals(mobj.getChannel())) {
            handleTransportChannel(mobj);
        } else if (mobj.getSocketID().equals(socket_id)) {
            Router.get().receiveMessage(socket_id,
                                        mobj.getChannel(),
                                        mobj.getMessage(),
                                        mobj.getFlags());
        }
        ***/    	
    }

    @Override
    public void onOpen(WebSocket.Connection new_socket) {
        socket = new_socket;
    }

    @Override
    public void onClose(int closeCode, String message) {
        wsobj_registry.remove(socket_id);
        socket_id = null;
    }

    private void handleTransportChannel(MessageObj mobj) {
        /* TODO: do proper UUID checks
         * see also: XHRServlet.handleTransportChannel()
         */
        if (mobj.getSocketID().equals("")) {
            /* Client wants socket ID
             * Because WSObj.create() already generated one,
             * there is nothing to do here.
             */
        }

        //Router.get().sendMessage(socket_id, Router.get().getInternalChannel(), "", "");
    }

    /* Factory stuff */
    static WSObj create() {
        String new_socket_id = java.util.UUID.randomUUID().toString();
        WSObj wsobj = new WSObj();
        wsobj.socket_id = new_socket_id;
        wsobj_registry.put(wsobj.socket_id, wsobj);
        return wsobj;
    }

    static boolean hasSocketID(String socket_id) {
        return wsobj_registry.containsKey(socket_id);
    }

    static void sendMessage(String socket_id, String message) {
        try {
            wsobj_registry.get(socket_id).socket.sendMessage(message);
        } catch (IOException e) {
            Logger.log(Logger.Level.WARN, "cannot send message", e);
        }
    }

    private String socket_id = null;
    private WebSocket.Connection socket = null;

    private static Map<String, WSObj> wsobj_registry;
    static {
        wsobj_registry = new HashMap<String, WSObj>();
    }
}

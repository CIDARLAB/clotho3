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
package org.clothocad.core.aspects.Router;

import java.io.PrintWriter;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.HashMap;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.json.JSONArray;
import org.json.JSONObject;

import org.clothocad.core.util.FileUtils;
import org.clothocad.core.util.Logger;

/**
 * @author Kelvin Li
 */
class XHRServlet extends HttpServlet {
    private XHRServlet() {
        send_buffers = new HashMap<String, LinkedList<JSONObject>>();
    }

    @Override
    public void doPost(HttpServletRequest request,
                       HttpServletResponse response)
                      throws ServletException, IOException {
        String received_data = getReceivedData(request);
        if (received_data == null) { return; }
        MessageObj mobj = MessageObj.fromJSONString(received_data);

        String new_socket_id = mobj.getSocketID();
        if (Router.get().getInternalChannel().equals(mobj.getChannel())) {
            new_socket_id = handleTransportChannel(mobj);
        } else if (hasSocketID(mobj.getSocketID())) {
            /* call parent */
            Router.get().receiveMessage(mobj.getSocketID(),
                                        mobj.getChannel(),
                                        mobj.getMessage(),
                                        mobj.getFlags());
        } /* otherwise drop the message */

        try {
            flushBuffer(new_socket_id, response);
        } catch (IOException e) {
            Logger.log(Logger.Level.WARN, "Cannot flush buffer", e);
        }
    }

    boolean hasSocketID(String socket_id) {
        return send_buffers.containsKey(socket_id);
    }

    void sendMessage(String socket_id, JSONObject message) {
        send_buffers.get(socket_id).add(message);
    }

    private String handleTransportChannel(MessageObj mobj) {
        /* TODO: polls come in on internal channel, 
         * but polls can have stale UUID's
         * need to do proper checking and re-issuing
         */
        String socket_id;
        if (mobj.getSocketID().equals("")) {
            socket_id = java.util.UUID.randomUUID().toString();
            send_buffers.put(socket_id, new LinkedList<JSONObject>());
            Router.get().sendMessage(socket_id, Router.get().getInternalChannel(), "", "");
        } else {
            socket_id = mobj.getSocketID();
        }

        return socket_id;
    }

    private void flushBuffer(String socket_id,
                             HttpServletResponse response)
                 throws IOException {
        List<JSONObject> buffer = send_buffers.get(socket_id);
        JSONArray output_json = new JSONArray(buffer);
        String output_data = output_json.toString();

        response.setStatus(HttpServletResponse.SC_OK);
        response.setContentType("text/plain");
        PrintWriter writer = response.getWriter();
        writer.println(output_data);
        writer.flush();
        writer.close();
        response.flushBuffer();
        buffer.clear();
    }

    private static String getReceivedData(HttpServletRequest request) {
        byte[] contents;
        try {
            contents = FileUtils.dumpInputStream(request.getInputStream());
        } catch (IOException e) {
            Logger.log(Logger.Level.WARN, "Cannot read incoming data", e);
            return null;
        }
        return new String(contents);
    }

    static XHRServlet get() {
        return singleton;
    }

    private static final XHRServlet singleton = new XHRServlet();

    private Map<String, LinkedList<JSONObject>> send_buffers;
}

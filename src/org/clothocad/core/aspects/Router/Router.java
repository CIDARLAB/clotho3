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

import java.util.Map;
import java.util.Set;
import javax.swing.JOptionPane;
import org.clothocad.core.aspects.Aspect;
import org.clothocad.core.aspects.Communicator.Communicator;
import org.clothocad.core.aspects.Communicator.Mind.Mind;
import org.clothocad.core.aspects.Communicator.Mind.PageMode;
import org.clothocad.core.aspects.Communicator.SendChannels;
import org.clothocad.core.aspects.Interpreter.AutoComplete;
import org.clothocad.core.aspects.Interpreter.Interpreter;
import org.clothocad.core.aspects.Logger;
import org.clothocad.core.datums.Doo;
import org.json.JSONArray;
import org.json.JSONObject;

/**
 * Low-level server<-->client communication layer
 *
 * Uses bidirectional WebSocket (WS) or XMLHttpRequest (XHR) connections
 * to pass messages. A message is a group of four ordered String
 * arguments: (socket_id, channel, message, flags).
 *
 * Serves web pages.
 *
 * Only Communicator should access Router.
 * 
 * Needs to match up with `web/scripts/clotholib/libsend.js`
 * 
 * @author Kelvin Li
 */
public class Router implements Aspect {
    public void sendMessage(String socket_id,
                            String channel,
                            String message,
                            String flags) {
        /* TODO: decide whether server-to-client flags are necessary */
        MessageObj mobj = MessageObj.fromComponents(socket_id,
                                                    channel,
                                                    message,
                                                    flags);
        if (WSObj.hasSocketID(socket_id)) {
            WSObj.sendMessage(socket_id, mobj.getJSON().toString());
        } else if (XHRServlet.get().hasSocketID(socket_id)) {
            XHRServlet.get().sendMessage(socket_id, mobj.getJSON());
        } else {
            Logger.log(Logger.Level.WARN, "Invalid socket ID");
        }
    }

    /* Called by XHRServlet and WSObj to
     * propagate incoming messages "upward"
     * 
     * JCA:  This is the first bottleneck entry point for a command coming from the Client
     * Previous to this call, the message could be handled by different mechnisms such as WebSocket
     * or Servelets, but they all get relayed here into this one spot, regardless of channels or 
     * anything like that.  All client traffic comes through this point.
     */
    void receiveMessage(String socket_id,
                        String channelStr,
                        String message,
                        String flags) {
        System.out.println("Router has received a message.");
        
        try {
            //
            //Instantiate the RouterDoo to manage this message
            //
            ClientChannel channel = ClientChannel.valueOf(channelStr);
            RouterDoo doo = new RouterDoo(socket_id, channel, message, flags);

            //
            //Fetch the Mind qassociated with this message
            //
            JSONObject json = new JSONObject(flags);
            String auth_key = json.getString("auth_key");
            /* TODO: validate auth_key */
            Mind mind = Communicator.get().getMind(auth_key);
            doo.mindId = mind.getId();
            
            //
            //Route the message via the channel enum
            //
            switch(channel) {
                case getTabList:
                    getTabList(doo, mind);
                    break;
                case linkPage:
                    linkPage(doo, mind);
                    break;
                case linkWidget:
                    linkWidget(doo, mind);
                    break;
                case serverEval:
                    serverEval(doo, mind);
                    break;
                case submitCommand:
                    submitCommand(doo, mind);
                    break;
                case unlinkPage:
                    unlinkPage(doo, mind);
                    break;
                case unlinkWidget:
                    unlinkWidget(doo, mind);
                    break;
                case updateQuery:
                    updateQuery(doo, mind);
                    break;
            }
        } catch(Exception err) {
            Logger.log(Logger.Level.WARN, "No handlers for " + channelStr);
        }
    }

    private void getTabList(RouterDoo doo, Mind mind) throws Exception {
        JSONArray out = new JSONArray();
        for (Map page_info : mind.getPageSummary()) {
            out.put(page_info);
        }
        Communicator.get().sendClientMessage(doo.socketId,
                                             SendChannels.showTabList,
                                             out.toString());
    }

    private void linkPage(RouterDoo doo, Mind mind) throws Exception {
            doo.message = new JSONObject(doo.messageStr);
            String ephemeral_id = doo.message.getString("ephemeral_link_page_id");
            String page_mode = doo.message.getString("page_mode");
Logger.log(Logger.Level.INFO, "linking w/ ephemeral ID: " + ephemeral_id);
            mind.linkPage(doo.socketId,
                          ephemeral_id,
//JCA:  this is ackward calls, should do it some other way
                          Enum.valueOf(PageMode.class, page_mode));
    }

    private void linkWidget(RouterDoo doo, Mind mind) throws Exception {
        Logger.log(Logger.Level.WARN, "linkWidget: Not implemented yet");
    }

    private void serverEval(RouterDoo doo, Mind mind) throws Exception {
        /**
         * Similar to "submitCommand" channel, but is intended to be used
         * by interactive widgets which need a simple serverEval function.
         * (e.g. to invoke assistants.)
         */
        if (!mind.eval(doo.socketId, doo.messageStr)) {
            /* TODO: tell client about this failure */
        }
    }
    
    /**
     * This is the method called when the user has clicked submit/search on the command bar
     * The message relayed here is an execution statement and interpreted as such, but it
     * needs to determine first whether it's native or computer text
     * 
     * @param doo
     * @param mind
     * @throws Exception 
     */
    private void submitCommand(RouterDoo doo, Mind mind) throws Exception {
        Logger.log(Logger.Level.INFO, "got message " + doo.messageStr);

        doo.message = new JSONObject(doo.messageStr);
        String cmd = doo.message.getString("command");

        //Try running the String; if it's formal code it will execute and return true.  If native returns false.
        if (!mind.runCommand(doo.socketId, cmd)) {
            disambiguate(mind, doo.socketId, cmd);
            completer.put(cmd);
        }
    }
    
        /**
         * If the command has failed to run already, then this interprets things as native
         * @param socket_id
         * @param cmd 
         */
        private void disambiguate(Mind mind, String socket_id, String nativeCmd) {
            Set<String> cmdResults = Interpreter.get().receiveNative(nativeCmd);

            if(cmdResults.isEmpty()) {
                /* TODO */
                //Communicator.get().say("(no search results available)");
                Logger.log(Logger.Level.WARN, "I got no search results");
            } else {
                //Send the search results

//TODO: TEMP ALERT THE USER
                Object[] possibilities = cmdResults.toArray();
                String assoc_cmd = (String)JOptionPane.showInputDialog(
                        null,
                        "Which to run?",
                        "Customized Dialog",
                        JOptionPane.PLAIN_MESSAGE,
                        null,
                        possibilities,
                        "ham");

                if (mind.runCommand(socket_id, assoc_cmd)) {
                    // learn to associate nativeCmd with assoc_cmd
                } else {
                    // tell the user that native cmd failed, so the overall action requested was aborted.
                    //Commicator.get().say("sorry man...no dice");
                }
                //IT SHOULD RETURN THE SEARCH RESULTS TO THE CLIENT, THEN LET THEM EXECUTE, BUT SINCE WE DON'T HAVE THAT NOW, IT WILL JUST RUN IT
                //I SEE A POTENTIAL ISSUE HERE OF NOT HAVING ACCESS TO THE MIND IN QUESTION
            }
        }

    private void unlinkPage(RouterDoo doo, Mind mind) throws Exception {
        mind.unlinkPage(doo.socketId);
    }

    private void unlinkWidget(RouterDoo doo, Mind mind) throws Exception {
        mind.unlinkWidget(doo.socketId, doo.messageStr);
    }
    
    /**
     * Handle query completions...so, this is invoked every time a key is pressed
     * to update the list of potential commands
     * @author jcanderson, Kelvin Li
     */
    private void updateQuery(RouterDoo doo, Mind mind) throws Exception {
        JSONObject out = getSuggestions(doo.messageStr);
        Communicator.get().sendClientMessage(
                                        doo.socketId,
                                        SendChannels.showQueryCompletions,
                                        out.toString());
    }
    
        private JSONObject getSuggestions(String message) throws Exception {
            JSONObject out = new JSONObject();
            JSONArray arr = wrapSuggestionsJSON(completer.getCompletions(message));
            out.put("query", message);
            out.put("suggestions", arr);
            return out;
        }

        private JSONArray wrapSuggestionsJSON(Iterable<String> suggestions) throws Exception {
            JSONArray arr = new JSONArray();
            for (String suggestion : suggestions) {
                JSONObject item = new JSONObject();
                item.put("command", suggestion);
                arr.put(item);
            }
            return arr;
        }
    
    /**
     * The Doo's that manage any Client-derived message.  Doo's handle even key
     * commands to avoid synchronization issues.
     */
    public class RouterDoo extends Doo {
        public RouterDoo(String sock, ClientChannel channel, String msg, String flags) {
            super(null, false);
            this.socketId = sock;
            this.channel = channel;
            this.messageStr = msg;
            this.flags = flags;
        }
        
        String socketId;
        ClientChannel channel;
        String messageStr;
        JSONObject message;
        String flags;
        String mindId;
    }

    String getInternalChannel() {
        return INTERNAL_CHANNEL;
    }

    private Router() {
        RouterServer.run();
    }
    
    public static enum ClientChannel {        
        getTabList,
        linkPage,
        linkWidget,
        serverEval,
        submitCommand,
        unlinkPage,
        unlinkWidget,
        updateQuery
    }

    public static Router get() {
        return singleton;
    }
    
    private final String INTERNAL_CHANNEL = "_transport";
    private final AutoComplete completer = new AutoComplete();

    private static Router singleton = new Router();
}

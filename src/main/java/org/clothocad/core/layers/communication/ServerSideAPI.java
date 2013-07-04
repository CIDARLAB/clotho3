    //clotho.show('CT-sample-view', '['static-admin-instance-is-uuid']', '{}'); clotho.show('CT-sample-view', '[]', '{}'); 
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

import com.fasterxml.jackson.core.JsonParseException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import javax.persistence.EntityNotFoundException;
import javax.persistence.NonUniqueResultException;
import javax.script.ScriptException;
import lombok.extern.slf4j.Slf4j;
import org.apache.shiro.authz.UnauthorizedException;
import org.bson.BSONObject;
import org.bson.types.ObjectId;
import org.clothocad.core.aspects.Interpreter.AutoComplete;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.aspects.Interpreter.Interpreter;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.Sharable;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.layers.communication.connection.ws.ClothoWebSocket;
import org.clothocad.core.layers.communication.mind.Mind;
import org.clothocad.core.layers.communication.mind.Widget;
import org.clothocad.core.util.FileUtils;
import org.clothocad.core.util.JSON;
import org.clothocad.model.Person;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * The ServerSideAPI relays the server methods that can be invoked by a client
 * on the server to the Aspects.
 *
 * It gets injected into scripting engines to handle conversion of language. It
 * goes into both the Mind's scripting engine and those of the Executor.
 *
 * It also is for securing the non-secure methods in the Aspects.
 *
 * So, it's just a wrapper really. It's all expressed in terms of Map<String,
 * Object>s and Lists to facilitate conversions. Each different context (Mind,
 * REST, websocket, java client, etc.) is expecting slightly different object
 * representations and synchronization models, so there is necessarily
 * interpretor logic in Router/ Communicator that handles this.
 *
 *
 * @author John Christopher Anderson
 */
@Slf4j
public final class ServerSideAPI {

    private static final Logger logger = LoggerFactory.getLogger(ServerSideAPI.class);
    private Persistor persistor;

    public ServerSideAPI(Persistor persistor) {
        this.persistor = persistor;
    }

// <editor-fold defaultstate="collapsed" desc="Human Interaction">      
    //JCA:  as of 6/6/2013 autcomplete works.  Wordlist is not persisted, but the completer does learn submitted phrases.
    public final void autocomplete(String userText) {
        List<String> completions = completer.getCompletions(userText);
        Message msg = new Message(Channel.autocomplete, completions);
        router.sendMessage(mind.getClientConnection(), msg);
    }

    //JCA:  works pushing a dummy message to the client, probably should be wrapped into get(...)
    public final String autocompleteDetail(String uuid) {
        try {
            Map<String, Object> msg = JSON.deserializeObject("{\"channel\":\"autocompleteDetail\",\"data\":{\"uuid\":\"1234567890\",\"text\":\"This is a command\",\"command\":\"clotho.run('230sdv-232', '18919e-18')\",\"versions\":[{\"uuid\":\"uuid123\",\"text\":\"Reverse Complement Tool\",\"author\":{\"uuid\":\"uuid_author_123\",\"name\":\"Joe Schmo\",\"email\":\"joe@schmo.com\",\"biography\":\"This is a biography about Joe Schmo. It's not too long. \"},\"description\":\"Aenean lacinia bibendum nulla sed consectetur. Cum sociis natoque penatibus et magnis dis parturient montes, nascetur ridiculus mus. Donec ullamcorper nulla non metus auctor fringilla. Maecenas faucibus mollis interdum. Etiam porta sem malesuada magna mollis euismod.\",\"usage\":{\"executed\":\"35\",\"successful\":\"27\",\"positive\":\"12\",\"negative\":\"3\"}},{\"uuid\":\"uuid456\",\"text\":\"pBca 1256\",\"author\":{\"uuid\":\"uuid_author_456\",\"name\":\"Chris Anderson\",\"email\":\"chris@anderson.com\",\"biography\":\"This is a biography about Chris Anderson. It's different than Joe's... It's a little longer. Yada yada yada. Here's some latin. It should get truncated on the server or we could write our own directive to handle truncating (easy). Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.\"},\"description\":\"Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.\",\"usage\":{\"executed\":\"8\",\"successful\":\"8\",\"positive\":\"6\",\"negative\":\"0\"}}]}}");
            return msg.get("data").toString();
        } catch (JsonParseException ex) {
            ex.printStackTrace();
            return null;
        }

    }

    //JCA:  as 0f 6/6/2013 submit seems to work
    public final void submit(String command) {
        //Resolve the arguments to a command string

        say(command, Severity.MUTED, null, true);

        if (!mind.runCommand(command)) {
            //disambiguate(command);  //JCA:  temporarily disabled for testing, also not fully hooked up
            return;
        }

        //If the command successfully executed, it gets retained
        mind.addLastCommand(Channel.submit, command);
    }

    public final void learn(Object data) {
        //might already be data?
        Map<String, Object> json = JSON.mappify(data);
        learn((String) json.get("userInput"), (Message) json.get("command"));
    }

    //clotho.learn("test1", "clotho.say('Hi there');")
    public final void learn(String userInput, Message command) {

        if (command == null) {
            command = mind.getLastCommands().get(0);
            if (command == null) {
                System.out.println("This is something that will happen if you haven't issued a previous command, it has nothing to learn");
                return;
            }
            Interpreter.get().learnNative(userInput, JSON.mappify(command));
        }
    }

    public final void login(String personRef, String password) {
    }

    public final void logout() {
    }

    public final boolean changePassword(String newPassword) {
        return true;
    }
// </editor-fold> 

    public final void clear() {
        mind.clear();
        say("The mind has been cleared", Severity.SUCCESS);
    }
// <editor-fold defaultstate="collapsed" desc="Logging and Messaging"> 
    //JCA:  as 0f 6/9/2013 say seems to work

    public final void say(Object obj) {
        if (obj instanceof String) {
            say((String) obj);
            return;
        }
        Map<String, Object> json = JSON.mappify(obj);
        Severity severity = json.containsKey("severity")
                ? Severity.valueOf(json.get("severity").toString())
                : Severity.NORMAL;
        String recipients = json.containsKey("recipients")
                ? json.get("recipients").toString()
                : null;
        boolean isUser = json.containsKey("isUser")
                ? Boolean.parseBoolean(json.get("isUser").toString())
                : false;
        say(json.get("message").toString(), severity, recipients, isUser);
    }

    public final void say(String message) {
        say(message, Severity.NORMAL, null, false);
    }

    /**
     *
     * @param message
     * @param severity "text-error", "text", "text-warning", "text-success" see
     * search-directives.js for 'from', is client or server
     */
    //JCA:  as 0f 6/9/2013 say seems to work
    public final void say(String message, Severity severity) {
//        System.out.println("say has : " + message);
        say(message, severity, null, false);

    }

    public final void say(String message, Severity severity, String recipients) {
//        System.out.println("say has : " + message);
        say(message, severity, recipients, false);

    }

    private final void say(String message, Severity severity, String recipients, boolean isUser) {

        //Resolve the recipients
        //XXX: doesn't currently handle multiple recipients
        //List<Sharable> listUsers = resolveToExistentSharablesList(recipients);
        String source;
        if (isUser) {
            source = "client";
        } else {
            source = "server";
        }
        Map<String, Object> data = new HashMap<>();
        data.put("text", message);
        data.put("from", source);
        data.put("class", severity);
        data.put("timestamp", new Date().getTime());

        Message msg = new Message(Channel.say, data);
        router.sendMessage(mind.getClientConnection(), msg);
    }

    //JCA:  Java side looks fine, but client code crashes browser
    //clotho.alert("this is an alert!");
    public final void alert(String message) {

        router.sendMessage(mind.getClientConnection(), new Message(Channel.alert, message));
    }

    //JCA:  This runs, and the message goes to the console.log spot.
    //clotho.log("I did some minipreps today");
    public final void log(String message) {
        log.debug("log has: {}", message);

        Message msg = new Message(Channel.log, message);

        router.sendMessage(mind.getClientConnection(), msg);
    }

    //Make note of this message in my notebook
    public final void note(String message) {
        System.out.println("I need to put this in your notebook, but i'm not implemented");
        //These have the same structure as a say, but they are stored.

        say("I've stored your note (but not really): " + message);
    }

// </editor-fold> 
// <editor-fold defaultstate="collapsed" desc="Data Manipulation"> 
    private void send(Message message) {
        router.sendMessage(mind.getClientConnection(), message);
    }

    public final void get(Object o, String requestId) {
        List<Map<String, Object>> returnData = new ArrayList<>();
        //list of selectors?
        if (o instanceof List) {
            for (Object obj : (List) o) {
                returnData.add(get(persistor.resolveSelector(obj.toString(), false)));
            }
        } else {
            returnData.add(get(persistor.resolveSelector(o.toString(), false)));
        }
        Message message = new Message(Channel.get, returnData, requestId);
        send(message);

    }

    //clotho.get("51b362f15076024ba019a642");  for a trail
    //clotho.get("51b29e7450765ced18af0d33");
    //JCA:  the object is requested, println'd, and collected in the clientside collector 6/8/2013
    //JCA:  result is also injected into the Mind's scriptengine as a proper object
    public final Map<String, Object> get(ObjectId id) {
        try {
            Map<String, Object> out = persistor.getAsJSON(id);
            say(String.format("Retrieved object #%s", id.toString()), Severity.SUCCESS);
            return out;

        } catch (UnauthorizedException e) {
            say(String.format("The current user does not have read access for object #%s", id.toString()), Severity.FAILURE);
            return null;
        } catch (EntityNotFoundException e) {
            say(String.format("No object with id %s", id.toString()), Severity.FAILURE);
            return null;
        }
    }

    private ObjectId resolveId(String id) {
        ObjectId uuid;
        try {
            uuid = new ObjectId(id);
        } catch (IllegalArgumentException e) {
            say(e.getMessage(), Severity.FAILURE);
            return null;
        }
        return uuid;
    }

    private void publish(Map obj) {
        router.publish(obj);
    }

    public final void set(Map<String, Object> values, String id) {
        try {

            if (values.get("id") == null) {
                say("No uuid provided", Severity.FAILURE);
                return;
            }

            ObjectId uuid = resolveId(values.get("id").toString());
            if (uuid == null) {
                return;
            }

            if (!persistor.has(uuid)) {
                say("No object with this id exists", Severity.FAILURE);
                return;
            }

            //Grab the object to be altered
            Map<String, Object> original = persistor.getAsJSON(uuid);

            //TODO: run validation check on saved values
            persistor.save(values);

            //Confirm that the new data is different than the old data
            Map<String, Object> altered = persistor.getAsJSON(uuid);
            if (original.toString().equals(altered.toString())) {
                say(String.format("Object #%s named %s was unmodified.", uuid.toString(), altered.get("name")), Severity.WARNING);
            }

            //Contact the user to notify them that they modified an object
            say(String.format("Successfully modified object #%s named %s", uuid.toString(), altered.get("name")), Severity.SUCCESS);

            //Relay the data change to listening clients
            publish(altered); //publish by uuid?
        } catch (UnauthorizedException e) {
            say(e.getMessage(), Severity.FAILURE);
        } catch (Exception e) {
            logAndSayError(String.format("Error setting %s: %s", values.toString(), e.getMessage()), e);
        }
    }

    public final void create(Object o, String requestId) {
        List<ObjectId> returnData = new ArrayList<>();
        //list of selectors?
        if (o instanceof List) {
            for (Object obj : (List) o) {
                returnData.add(create(JSON.mappify(obj)));
            }
        } else {
            returnData.add(create(JSON.mappify(o)));
        }
        Message message = new Message(Channel.create, returnData, requestId);
        send(message);
    }

    public final ObjectId create(Map<String, Object> obj) {

        try {
            //Confirm that there is no pre-existing object with this uuid
            String idKey = null;
            if (obj.containsKey("id")) idKey = "id";
            if (obj.containsKey("_id")) idKey = "_id";
            
            
            if (idKey != null) {
                
                ObjectId uuid = resolveId(obj.get(idKey).toString());
                if (uuid == null) {
                    return null;
                }
                if (persistor.has(uuid)) {
                    say("An object with the uuid " + uuid + " already exists.  No object was created.", Severity.FAILURE);
                    return null;
                }
                
                obj.put(idKey, new ObjectId(obj.get(idKey).toString()));
            }

            ObjectId id = persistor.save(obj);


            //TODO: Relay the data change to listening clients
            System.out.println("Ernst, this needs to be implemented here too.  Push object via pubsub.");

            //Return the JSON of the new object as a String
            say(String.format("Created object #%s named %s", id.toString(), obj.get("name")), Severity.SUCCESS);
            return id;
        } catch (UnauthorizedException e) {
            say("The current user does not have write access for this domain", Severity.FAILURE);
            return null;
        } 
        
        //catch (Exception e) {
        //    logAndSayError(String.format("Error creating %s: %s", obj.toString(), e.getMessage()), e);
        //    return null;
        //}
    }

    private void logAndSayError(String message, Exception e) {
        log.error(message, e);
        say(message, Severity.FAILURE);
    }

    public final void destroy(Object o, String requestId) {
        //list of selectors?
        if (o instanceof List) {
            for (Object obj : (List) o) {
                destroy(persistor.resolveSelector(obj.toString(), false));
            }
        } else {
            destroy(persistor.resolveSelector(o.toString(), false));
        }
    }

    //JCA:  as of 6/9/2013 works
    public final void destroy(ObjectId id) {
        try {
            try {
                persistor.delete(id);
            } catch (UnauthorizedException e) {
                say(e.getMessage(), Severity.FAILURE);
            }
            say(String.format("Destroyed object #%s", id.toString()), Severity.SUCCESS);
        } catch (Exception e) {
            logAndSayError(String.format("Error destroying %s: %s", id.toString(), e.getMessage()), e);
        }
    }

    public final void query(Map<String, Object> spec, String requestId) {
        List<Map<String, Object>> objs;

        //XXX: demo hack to resolve schema smartly
        if (spec.containsKey("schema")) {
            //figure out what the schema actually is
            try {
                Map<String, Object> schema = persistor.getAsJSON(persistor.resolveSelector(spec.get("schema").toString(), false));
                String schemaName = schema.get("binaryName").toString(); //try and fallback to name name?
                spec.remove("schema");
                spec.put("className", schemaName);
            } catch (EntityNotFoundException e) {
                logAndSayError(String.format("No schema found for selector %s", spec.get("schema").toString()), e);
                objs = new ArrayList<>();
                Message msg = new Message(Channel.query, objs, requestId);
                router.sendMessage(mind.getClientConnection(), msg);
                return;
            }
        }

        try {
            //Relay the query to Persistor and return the hits
            objs = persistor.findAsBSON(spec);
            say("Found " + objs.size() + " objects that satisfy your query", Severity.SUCCESS);
            Message msg = new Message(Channel.query, objs, requestId);
            router.sendMessage(mind.getClientConnection(), msg);
        } catch (Exception e) {
            logAndSayError(String.format("Error querying %s: %s", spec.toString(), e.getMessage()), e);
        }

    }
// </editor-fold> 

// <editor-fold defaultstate="collapsed" desc="Execution"> 
    public final void run(Object o, String id) throws ScriptException {
        Map<String, Object> data = JSON.mappify(o);

        Function function = persistor.get(Function.class, persistor.resolveSelector(data.get("function").toString(), Function.class, false));
        Map<String, Object> arguments = data.containsKey("arguments")
                ? JSON.mappify(data.get("arguments"))
                : null;
        Object result = run(function, arguments);
        if (!result.equals(Function.NoResult.class)) {
            //TODO: Map<String, Object> reply - we need to have a way to designate which request we are responding to
            Message message = new Message(Channel.run, result, id);
            send(message);
        }

    }

    public final Object run(Function function, Map<String, Object> args) throws ScriptException {
        return function.execute(args);
    }

    /**
     * Relay method for receiving a show command
     *
     * @param viewRef
     * @param sharables
     * @param sloppy_args
     */
    public final void show(String viewIod, String sharables, String position, String recipients) {
        System.out.println("Show is called " + position);
        /**
         * *
         * try {
         *
         * //IF THE PROCESS WAS TIED TO A COMMAND-INITIATED PROCESS THERE WOULD
         * BE A PARENT DOO THAT NEEDS TO BE FOUND. Doo parentDoo =
         * Hopper.get().extract(null); ShowDoo doo = new ShowDoo(parentDoo);
         *
         * //Gather up referenced objects View view = (View)
         * resolveToObjBase(viewId); List<Sharable> shareList =
         * resolveToSharables(sharables);
         *
         * //Create the widget and put it on its page Page targetPage =
         * mind.getConfig().getPage(socket_id); Widget widget = new
         * Widget(targetPage, view); targetPage.addWidget(widget);
         *
         * //Just for record keeping doo.viewId = view.getId(); doo.widgetId =
         * widget.getId(); doo.collectSharables(shareList);
         *
         * //Create all the commands for the client doo.commandMessageArray =
         * new ArrayList(); doo.commandMessageArray.put(makeCollect(shareList));
         * doo.commandMessageArray.put(makeShowWidget(widget, position));
         * doo.commandMessageArray.put(makeUpdate(widget, shareList,
         * socket_id)); doo.commandMessageArray.put(makeCallback(doo.getId(),
         * socket_id));
         *
         * //Put the doo into the hopper to await a callback
         * Hopper.get().add(doo);
         *
         * //Send the commands Communicator.get().sendClientMessage(socket_id,
         * SendChannels.commandList, doo.commandMessageArray.toString());
         *
         * //Save everything whose state was changed //JCA: THIS NEEDS A
         * CALLBACK/FAILURE RESPONSE THAT REVERTS THIS (EVENTUALLY)
         * Persistor.get().persistObjBase(mind); } catch (Exception e) {
         * Logger.log(Logger.Level.WARN, "", e); e.printStackTrace(); //REVERT
         * THE MIND //SEND CLIENT MESSAGE TELLING THAT IT FAILED } *
         */
    }

    public final void startTrail(String trailRef) {
        /**
         * *
         * try { //Grab a Doo if something is awaiting one Doo parentDoo =
         * Hopper.get().extract(null); Trail trail = (Trail)
         * resolveToObjBase(trailRef); Person student = this.getPerson();
         *
         * Proctor.get().initiateTrail(student, trail, parentDoo);
         *
         * } catch (Exception e) { Logger.log(Logger.Level.WARN, "", e);
         * e.printStackTrace(); } *
         */
    }

    /**
     * Invoke the hard-coded editor appropriate for this particular object *For
     * an Instance, invoke the view set as the default view for this type of
     * data. That's one option. The other is to build up an editor
     * programmatically from the components. Not sure. Actually, I like this
     * second option best. Maybe. I don't know. Second option is more secure.
     *
     * Both options need to be available, I think.
     *
     * @param sharableRef
     */
    public final void edit(String sharableRef) {
        /**
         * try { Sharable sharable = (Sharable) resolveToObjBase(sharableRef);
         *
         * switch (sharable.type()) { case SCHEMA: Schema schema = (Schema)
         * sharable; //Pop up the page, push in data return; case INSTANCE:
         * return; default: return; } //Pop up a new window of } catch
         * (Exception err) { err.printStackTrace(); }
         *
         */
    }

    public final void listen(String args) {
        say("not yet implemented", Severity.FAILURE);
    }

    public final void unlisten(String data) {
        say("not yet implemented", Severity.FAILURE);
    }

// </editor-fold> 
    // <editor-fold defaultstate="collapsed" desc="Setup"> 
    public ServerSideAPI(Mind mind, Persistor persistor) {
        this.mind = mind;
        this.persistor = persistor;
    }

    public void setMind(Mind mind) {
        this.mind = mind;
    }
// </editor-fold>

//    public final void test() {
//        try {
//            System.out.println("Test has been invoked");
//            Map<String, Object> trail = new HashMap<String, Object>("{\"uuid\":\"trail_123\",\"title\":\"Biosafety Module\",\"author\":\"UC Berkeley\",\"description\":\"<blockquote><p>This is a module on Biosafety. You'll learn about Cras sit amet nibh libero, in gravida nulla. Nulla vel metus scelerisque ante sollicitudin commodo. Cras purus odio, vestibulum in vulputate at, tempus viverra turpis.</p><p><small>Maxwell Bates</small></p></blockquote>\",\"contents\":[{\"module_title\":\"The Basics\",\"pavers\":[{\"paver_title\":\"Introduction\",\"type\":\"template\",\"template\":\"/app/partials/trail_123_1.html\"}]},{\"module_title\":\"Biosafety Levels\",\"pavers\":[{\"paver_title\":\"Introduction\"},{\"paver_title\":\"Biosafety Level 1-2\"},{\"paver_title\":\"Biosafety Level 3-4\"}]},{\"module_title\":\"Assessment\",\"pavers\":[{\"paver_title\":\"Review\"}],\"assessment\":[{\"type\":\"quiz\",\"title\":\"Final Quiz\"}]}]}");
//            Trail i = new Trail("Biosafety Module", trail);
//
//            System.out.println("Trail i has been created: " + i.toString());
//            
//            Persistor.get().save(i);
//            ObjectId id = i.getUUID();
//            
//            String uuid = id.toString();
//
//            say(id.toString());
//            
//            ObjBase reclaimed =  Collector.get().getObjBase(uuid);
//            
//            System.out.println("After re-retrieval from db I have: " + reclaimed.toString());
//        } catch (JSONException ex) {
//            java.util.logging.Logger.getLogger(ServerSideAPI.class.getName()).log(Level.SEVERE, null, ex);
//        }
//    }
    /**
     * Return the Person object associated with the Mind into which this SS API
     * has been injected
     *
     * @return
     */
    private Person getPerson() {
        return mind.getPerson();
        /**
         * if (person != null) { return person; } return
         * Person.getById(mind.getPersonId());
         *
         */
    }

    public final void newPage() throws Exception {
        /**
         * *
         * try { List commandMessageArray = new ArrayList();
         * commandMessageArray.put(makeNewPage(UUID.randomUUID().toString()));
         *
         * // // // //Put the doo into the hopper to await a callback //
         * Hopper.get().add(doo);
         *
         * //Send the commands Communicator.get().sendClientMessage(socket_id,
         * SendChannels.commandList, commandMessageArray.toString());
         *
         * //Save everything whose state was changed //JCA: THIS NEEDS A
         * CALLBACK/FAILURE RESPONSE THAT REVERTS THIS (EVENTUALLY)
         * Persistor.get().persistObjBase(mind); } catch (Exception e) {
         * Logger.log(Logger.Level.WARN, "", e); e.printStackTrace(); //REVERT
         * THE MIND //SEND CLIENT MESSAGE TELLING THAT IT FAILED } *
         */
    }

    /**
     * //if a command has a callback, handle it here //minimally removes the
     * doo from the dooHopper //jca: i don't know that i like this being in the
     * api like this. it really //should be on a more secure channel because
     * this all deals with Doos, which //need to be accurately tracked. if you
     * can directly call this callback from //the api, you can just type it into
     * the commmand bar and lie to clotho, which //is no good. so, yeah, it
     * shouldn't be implemented like this, but we can fix that //in the future
     * when we pick through the client-->server communications standard. //also,
     * this is a mind-agnostic thing -- this is what gets called just to confirm
     * //that a process has completed, and everything should probably call this.
     * The doo //hopper maybe should even be an Aspect with this being relayed
     * there. Yeah. That's //it.
     */
    public final void notify(String dooID) {
        /**
         * System.out.println("calling back!"); //Grab the Doo and perhaps
         * perform any commands embedded in there
         *
         * //Terminate it from the Hopper Hopper.get().terminate(dooID);
         *
         */
    }

    // <editor-fold defaultstate="collapsed" desc="Command Factories"> 
    /**
     * Instantiate a new workspace page. Only the Proctor can construct
     * arguments for TRAILS or specific editors. Security thing.
     *
     * @return
     * @throws Exception
     */
    public static Map<String, Object> makeNewPage(String page_id) throws Exception {
        Map<String, Object> out = new HashMap<>();
        out.put("mode", "WORKSPACE");
        out.put("ephemeral_link_page_id", page_id);
        out.put("command", "addPage");
        return out;
    }

    /**
     * Instruct the client to display the GUI of a view and position it
     * according to positioning parameters
     *
     * @param viewId
     * @param sloppy_args
     * @return
     * @throws Exception
     */
    public Map<String, Object> makeShowWidget(Widget widget, String position) throws Exception {
//        Map<String, Object> out = new HashMap<String, Object>();
//        Map<String, Object> positioninfo = resolveToMap<String, Object>(position);
//        View view = widget.getView();
//
//        //Load the html and js scripts into a Map<String, Object>, insert the uuid for the widget
//        String widgetId = widget.getId();
//        System.out.println("graphics script: "+view.getGraphicsScript());
//        String html = replaceWidgetId(view.getGraphicsScript(), widgetId);
//        String onshow = replaceWidgetId(view.getOnShowScript(), widgetId);
//
//        out.put("widget_id", widgetId);
//        out.put("content", html);
//        out.put("on_show", onshow);
//        out.put("parent_widget_id", "widget_space_root");
//        out.put("command", "showWidget");
//
//        return out;
        return null;
    }

    /**
     * Update instructs the client to call update on a widget The data that goes
     * into the widget, however, is sent in a Collect command
     *
     * @param view
     * @return
     * @throws Exception
     */
    public static Map<String, Object> makeUpdate(Widget widget, List<Sharable> sharables, String socket_id) throws Exception {
        Map<String, Object> updateCmd = new HashMap<>();
        //Createhte token:uuid map that associates the piece of data with its token used by the widget
        List<ClothoField> fields = widget.getView().getInputArguments();
        Map<String, Object> inputs = new HashMap<>();
        for (int i = 0; i < fields.size(); i++) {
            ClothoField field = fields.get(i);
            inputs.put(field.getName(), sharables.get(i).getUUID());
        }
        updateCmd.put("inputArgs", inputs);
        updateCmd.put("script", replaceWidgetId(widget.getView().getOnUpdateScript(), "")); //how do i get widgetID's again?
        updateCmd.put("command", "update");
//            
        return updateCmd;
    }

    public static Map<String, Object> makeCallback(String dooId, String socket_id) throws Exception {
        Map<String, Object> callBackCmd = new HashMap<>();
        callBackCmd.put("socketID", socket_id);
        callBackCmd.put("dooID", dooId);
        callBackCmd.put("command", "callback");

        return callBackCmd;
    }
// </editor-fold> 

    // <editor-fold defaultstate="collapsed" desc="Utility Methods"> 
    /**
     * Replace the _widget_id phrases with the actual uuid
     *
     * @param script
     * @param widgetIdPrefix
     * @return
     */
    public static final String replaceWidgetId(String script, String widgetIdPrefix) {
        try {
            return XMLParser.addPrefixToTagAttribute(script, "id", widgetIdPrefix);
        } catch (Exception ex) {
            logger.error("", ex);
        }
        return null;
    }

    // </editor-fold> 
    /**
     * ****** VARIABLES *******
     */
    //private transient Person person;
    private static final AutoComplete completer = new AutoComplete();
    private Mind mind;
    private static Router router;

    static {
        router = Router.get();
    }

    public enum Severity {

        SUCCESS,
        WARNING,
        FAILURE,
        NORMAL,
        MUTED
    }
}

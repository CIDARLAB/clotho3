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
package org.clothocad.core.communication;

import org.clothocad.core.util.XMLParser;
import com.fasterxml.jackson.core.JsonParseException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import javax.persistence.EntityNotFoundException;
import javax.script.ScriptException;
import lombok.extern.slf4j.Slf4j;
import org.apache.shiro.SecurityUtils;
import org.apache.shiro.authc.AuthenticationException;
import org.apache.shiro.authc.UsernamePasswordToken;
import org.apache.shiro.authz.UnauthorizedException;
import org.apache.shiro.subject.Subject;
import org.bson.types.ObjectId;
import org.clothocad.core.aspects.Interpreter.AutoComplete;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.aspects.Interpreter.Interpreter;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.Module;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.execution.Mind;
import org.clothocad.core.communication.mind.Widget;
import org.clothocad.core.schema.BuiltInSchema;
import org.clothocad.core.schema.ReflectionUtils;
import org.clothocad.core.util.JSON;
import org.clothocad.model.Person;

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
public class ServerSideAPI {

    private static final AutoComplete completer = new AutoComplete();
    private final Router router;
    private final Persistor persistor;
    private final String requestId;
    private final Mind mind;

    public ServerSideAPI(Mind mind, Persistor persistor, Router router, String requestId) {
        this.persistor = persistor;
        this.mind = mind;
        this.requestId = requestId;
        this.router = router;
    }

// <editor-fold defaultstate="collapsed" desc="Human Interaction">      
    //JCA:  as of 6/6/2013 autcomplete works.  Wordlist is not persisted, but the completer does learn submitted phrases.
    public final void autocomplete(String userText) {
        List<String> completions = completer.getCompletions(userText);
        Message msg = new Message(Channel.autocomplete, completions);
        router.sendMessage(mind.getConnection(), msg);
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
    public final Object submit(String command) {
        //Resolve the arguments to a command string
        //say(command, Severity.MUTED, null, true);
        try {
            Object returnValue = mind.runCommand(command, new ScriptAPI(mind, persistor, router, requestId));
            //If the command successfully executed, it gets retained
            mind.addLastCommand(Channel.submit, command);
            return returnValue;
        } catch (ScriptException ex) {
            //disambiguate(command);  //JCA:  temporarily disabled for testing, also not fully hooked up
            logAndSayError("Error while executing script: " + ex.getMessage(), ex);
            return Void.TYPE;
        }
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

    public final boolean login(String username, String password) {
        try {
            SecurityUtils.getSubject().login(new UsernamePasswordToken(username, password));
            say("Welcome, " + username, Severity.SUCCESS);
            log.info("User {} logged in", username);
            return true;
        } catch (AuthenticationException e) {
            logAndSayError("Authentication attempt failed for username " + username, e);
            return false;
        }
    }

    public final boolean logout() {
        if (SecurityUtils.getSubject().isAuthenticated()) {
            String username = SecurityUtils.getSubject().getPrincipal().toString();
            mind.setUsername(username);
            persistor.save(mind);
            SecurityUtils.getSubject().logout();
            say("Logged out", Severity.SUCCESS);
            log.info("User {} logged out", username);
            return true;
        }

        say("You are not logged in", Severity.WARNING);
        return false;
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

    protected void say(String message, Severity severity, String recipients, boolean isUser) {

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

        Message msg = new Message(Channel.say, data, requestId);
        router.sendMessage(mind.getConnection(), msg);
    }

    //JCA:  Java side looks fine, but client code crashes browser
    //clotho.alert("this is an alert!");
    public final void alert(String message) {

        router.sendMessage(mind.getConnection(), new Message(Channel.alert, message));
    }

    //JCA:  This runs, and the message goes to the console.log spot.
    //clotho.log("I did some minipreps today");
    public final void log(String message) {
        log.debug("log has: {}", message);

        Message msg = new Message(Channel.log, message);

        router.sendMessage(mind.getConnection(), msg);
    }

    //Make note of this message in my notebook
    public final void note(String message) {
        System.out.println("I need to put this in your notebook, but i'm not implemented");
        //These have the same structure as a say, but they are stored.

        say("I've stored your note (but not really): " + message);
    }

// </editor-fold> 
// <editor-fold defaultstate="collapsed" desc="Data Manipulation"> 
    protected void send(Message message) {
        router.sendMessage(mind.getConnection(), message);
    }

    private Object unwrap(Object o) {
        while (o instanceof Collection) {
            Iterator iterator = ((Collection) o).iterator();
            if (iterator.hasNext()) {
                o = iterator.next();
            } else {
                throw new IllegalArgumentException("Cannot unwrap an empty collection");
            }
        }

        return o;
    }

    public Map<String, Object> get(Object o) {
        o = unwrap(o);
        return get(persistor.resolveSelector(o.toString(), false));
    }

    public Map<String, Object> get(ObjectId id) {
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

    public final List<Map<String, Object>> getAll(List objects) {
        List<Map<String, Object>> returnData = new ArrayList<>();
        for (Object obj : objects) {
            returnData.add(get(persistor.resolveSelector(obj.toString(), false)));
        }
        return returnData;
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

    public final List<ObjectId> setAll(List<Map<String, Object>> values) {
        List<ObjectId> out = new ArrayList<>();
        for (Map<String, Object> spec : values) {
            out.add(set(spec));
        }
        return out;
    }

    public final ObjectId set(Map<String, Object> values) {
        try {

            if (values.get("id") == null) {
                say("set: No uuid provided", Severity.WARNING);
                return create(values);
            }

            ObjectId uuid = resolveId(values.get("id").toString());
            if (uuid == null) {
                return null;
            }

            if (!persistor.has(uuid)) {
                say("set: No object with this id exists", Severity.WARNING);
                return create(values);
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
            return uuid;
        } catch (UnauthorizedException e) {
            say(e.getMessage(), Severity.FAILURE);
            return null;
        } catch (Exception e) {
            logAndSayError(String.format("Error setting %s: %s", values.toString(), e.getMessage()), e);
            return null;
        }
    }

    public final String create(Object o) {
        return create(JSON.mappify(o)).toString();
    }

    //TODO: some global solution for jsonifying ObjectIds
    public final List<String> createAll(List<Object> objects) {
        List<String> returnData = new ArrayList<>();
        //list of selectors?
        for (Object obj : objects) {
            returnData.add(create(JSON.mappify(obj)).toString());
        }
        return returnData;
    }

    public ObjectId create(Map<String, Object> obj) {

        try {
            //Confirm that there is no pre-existing object with this uuid
            String idKey = null;
            if (obj.containsKey("id")) {
                idKey = "id";
            }
            if (obj.containsKey("_id")) {
                idKey = "_id";
            }


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
            //TODO: create sets author to current user
            ObjectId id = persistor.save(obj);

            //TODO: create as java object and save to force validation and trigger post-create methods

            try {
                ObjBase realObject = persistor.get(ObjBase.class, id);
                persistor.save(realObject);
            } catch (Exception e) {
                logAndSayError("could not validate object - possibly malformed data?", e);
            }



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

    protected void logAndSayError(String message, Exception e) {
        log.error(message, e);
        say(message, Severity.FAILURE);
    }

    public final ObjectId destroy(Object id) {
        ObjectId resolvedId = persistor.resolveSelector(id.toString(), false);
        if (resolvedId == null) {
            return null;
        }
        try {
            try {
                persistor.delete(resolvedId);
            } catch (UnauthorizedException e) {
                say(e.getMessage(), Severity.FAILURE);
            }
            say(String.format("Destroyed object #%s", resolvedId.toString()), Severity.SUCCESS);
            return resolvedId;
        } catch (Exception e) {
            logAndSayError(String.format("Error destroying %s: %s", id.toString(), e.getMessage()), e);
            return null;
        }
    }

    public final List<ObjectId> destroyAll(List<Object> objects) {
        List<ObjectId> out = new ArrayList<>();
        for (Object obj : objects) {
            out.add(destroy(obj));
        }
        return out;
    }

    public List<Map<String, Object>> query(Map<String, Object> spec) {
        List<Map<String, Object>> objs;
        //in case we are handed an immutable map, make a copy
        spec = new HashMap(spec);

        //XXX: demo hack to resolve schema smartly
        if (spec.containsKey("schema")) {
            //figure out what the schema actually is
            try {
                Map<String, Object> schema = persistor.getAsJSON(persistor.resolveSelector(spec.get("schema").toString(), false));
                String binaryName = (String) schema.get("binaryName"); //try and fallback to name name?
                //if (schemaName == null || ){ 
                //    schemaName = schema.get("name").toString();
                //}
                String textName = (String) schema.get("name");
                
                if (!(schema instanceof BuiltInSchema)) textName = (String) schema.get("name");
                
                spec.remove("schema");
                if (binaryName != null && textName != null){
                    Map<String,Object> query = new HashMap();
                    List<String> names = new ArrayList<>();
                    names.add(textName);
                    names.add(binaryName);
                    query.put("$in", names);
                    spec.put("className", query);
                } else if (binaryName != null){
                    spec.put("className", binaryName);
                } else if (textName != null){
                    spec.put("className", textName);
                } else {
                    spec.put("className", ((ObjBase) schema).getUUID());
                }
            } catch (EntityNotFoundException e) {
                //maybe already full name?
                logAndSayError(String.format("No schema found for selector %s", spec.get("schema").toString()), e);
                objs = new ArrayList<>();
                return objs;
            }
        }

        try {
            //Relay the query to Persistor and return the hits
            objs = persistor.findAsBSON(spec);
            say("Found " + objs.size() + " objects that satisfy your query", Severity.SUCCESS);
            return objs;
        } catch (Exception e) {
            logAndSayError(String.format("Error querying %s: %s", spec.toString(), e.getMessage()), e);
            e.printStackTrace();
            return new ArrayList<>();
        }

    }
// </editor-fold> 

// <editor-fold defaultstate="collapsed" desc="Execution"> 
    //TODO: needs serious cleaning up
    public final Object run(Object o) throws ScriptException, IllegalAccessException, IllegalArgumentException, InvocationTargetException {
        Map<String, Object> data = JSON.mappify(o);
        List<Object> args;

        try {
            args = (List) data.get("args");
        } catch (ClassCastException e) {
            say("Arguments must be a list or array");
            return null;
        }

        if (data.containsKey("id")) {
            //XXX:(ugh ugh) end-run if *Function
            Map<String, Object> functionData = persistor.getAsJSON(persistor.resolveSelector(data.get("id").toString(), true));
            if (functionData.containsKey("schema") && functionData.get("schema").toString().endsWith("Function")) {
                try {
                    Function function = persistor.get(Function.class, persistor.resolveSelector(data.get("id").toString(), true));

                    return mind.invoke(function, args, new ScriptAPI(mind, persistor, router, requestId));
                } catch (ScriptException e) {
                    logAndSayError("Script Exception thrown: " + e.getMessage(), e);
                    return Void.TYPE;
                } catch (NoSuchMethodException ex) {
                    logAndSayError("No such function found", ex);
                    return Void.TYPE;
                }
            }
            //XXX: this whole function is still a mess
            if (functionData.containsKey("schema") && functionData.get("schema").toString().endsWith("Module")) {
                try {
                    Module module = persistor.get(Module.class, persistor.resolveSelector(data.get("id").toString(), true));

                    return mind.invokeMethod(module, data.get("function").toString(), args, new ScriptAPI(mind, persistor, router, requestId));
                } catch (ScriptException e) {
                    logAndSayError("Script Exception thrown: " + e.getMessage(), e);
                    return Void.TYPE;
                } catch (NoSuchMethodException ex) {
                    logAndSayError("No such method found", ex);
                    return Void.TYPE;
                }
            }

            for (int i = 0; i < args.size(); i++) {
                try {
                    ObjectId id = new ObjectId(args.get(i).toString());
                    args.set(i, persistor.get(ObjBase.class, id));
                } catch (IllegalArgumentException e) {
                }
            }
            //reflectively (ugh) run function of instance
            ObjBase instance = persistor.get(ObjBase.class, new ObjectId(data.get("id").toString()));



            Method method = ReflectionUtils.findMethodNamed(data.get("function").toString(), args.size(), instance.getClass());
            Object result = method.invoke(instance, args.toArray());
            if (method.getReturnType().equals(Void.TYPE)) {
                return Void.TYPE;
            }
            List<Object> results = new ArrayList<>();
            if (result instanceof Iterable) {
                for (Object r : ((Iterable) result)) {
                    if (r instanceof ObjBase) {
                        results.add(persistor.save((ObjBase) r));
                    } else {
                        results.add(r);
                    }
                }
            } else {
                if (result instanceof ObjBase) {
                    results.add(persistor.save((ObjBase) result));
                } else {
                    results.add(result);
                }
            }
            Message message = new Message(Channel.run, persistor.toJSON(results), requestId);
            send(message);
            return Void.TYPE;
        }

        Function function = persistor.get(Function.class, persistor.resolveSelector(data.get("function").toString(), Function.class, false));
        List arguments;
        try {
            arguments = data.containsKey("arguments")
                    ? JSON.deserializeList(data.get("arguments").toString())
                    : null;
        } catch (JsonParseException ex) {
            logAndSayError("malformed arguments", ex);
            return Void.TYPE;
        }
        Object result = run(function, arguments);
        if (!result.equals(Void.TYPE)) {
            Message message = new Message(Channel.run, result, requestId);
            send(message);
        }

        return Void.TYPE;

    }

    public final Object run(Function function, List<Object> args) throws ScriptException {


        return mind.evalFunction(function.getCode(), function.getName(), args, new ScriptAPI(mind, persistor, router, requestId));
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
            log.error("", ex);
        }
        return null;
    }

    Object queryOne(Map<String, Object> query) {
        List result = query(query);
        if (result.isEmpty()) {
            throw new EntityNotFoundException();
        }
        return result.get(0);
    }

    public enum Severity {

        SUCCESS,
        WARNING,
        FAILURE,
        NORMAL,
        MUTED
    }
}
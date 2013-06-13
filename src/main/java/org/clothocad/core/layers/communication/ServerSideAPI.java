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

import com.mongodb.BasicDBObject;
import com.mongodb.DBObject;
import com.mongodb.util.JSON;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import javax.script.ScriptEngine;
import javax.swing.JOptionPane;
import org.bson.BSONObject;
import org.bson.types.ObjectId;
import org.clothocad.core.aspects.Authenticator;
import org.clothocad.core.aspects.Collector;
import org.clothocad.core.aspects.Interpreter.AutoComplete;
import org.clothocad.core.aspects.Logger;
import org.clothocad.core.aspects.Persistor;
import org.clothocad.core.aspects.Interpreter.Interpreter;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.Sharable;
import org.clothocad.core.datums.View;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.layers.communication.connection.ws.ClothoWebSocket;
import org.clothocad.core.layers.communication.mind.Mind;
import org.clothocad.core.layers.communication.mind.Widget;
import org.clothocad.core.layers.persistence.mongodb.MongoDBConnection;
import org.clothocad.model.Person;
import org.clothocad.model.Trail;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;
import org.json.JSONTokener;

/**
 * The ServerSideAPI relays the server methods that can be invoked by a client on
 * the server to the Aspects.
 * 
 * It gets injected into scripting engines to handle conversion of language.  It
 * goes into both the Mind's scripting engine and those of the Executor.
 * 
 * It also is for securing the non-secure methods in the Aspects.
 * 
 * So, it's just a wrapper really.  It's all expressed in terms of JSONObjects
 * and JSONArrays to facilitate conversions.  Each different context (Mind, REST, 
 * websocket, java client, etc.) is expecting slightly different object representations
 * and synchronization models, so there is necessarily interpretor logic in Router/
 * Communicator that handles this.
 * 
 * I started implementing this at the level of Sharables, but that was super hard.  I
 * then went to doing it with JSONObject and JSONArray and got stuck on implementing 
 * create.  In the end, everything goes through a String representation somewhere, so
 * that will be the interface -- String representations of JSON
 * 
 * @author John Christopher Anderson
 */
public final class ServerSideAPI {

// <editor-fold defaultstate="collapsed" desc="Human Interaction">      
    //JCA:  as of 6/6/2013 autcomplete works.  Wordlist is not persisted, but the completer does learn submitted phrases.
    public final String autocomplete(String userText) {
        try {
            ArrayList<String> completions = completer.getCompletions(userText);
            JSONArray data = new JSONArray(); //The list of JSONObjects each for an autocomplete
            for(String str : completions) {
                JSONObject obj = new JSONObject();
                obj.put("text", str);
                obj.put("uuid", "irrelevent_uuid");
                obj.put("type", "phrase");
                data.put(obj);
            }
            
            JSONObject msg = new JSONObject();
            msg.put("channel", "autocomplete");
            msg.put("data", data);
            Router.get().sendMessage(mind.getClientConnection(), msg);
            return data.toString();
        } catch(Exception err) {
            err.printStackTrace();
            return null;
        }
    }
    
    //JCA:  works pushing a dummy message to the client, probably should be wrapped into get(...)
    public final String autocompleteDetail(String uuid) {
        try {
            JSONObject msg = new JSONObject("{\"channel\":\"autocompleteDetail\",\"data\":{\"uuid\":\"1234567890\",\"text\":\"This is a command\",\"command\":\"clotho.run('230sdv-232', '18919e-18')\",\"versions\":[{\"uuid\":\"uuid123\",\"text\":\"Reverse Complement Tool\",\"author\":{\"uuid\":\"uuid_author_123\",\"name\":\"Joe Schmo\",\"email\":\"joe@schmo.com\",\"biography\":\"This is a biography about Joe Schmo. It's not too long. \"},\"description\":\"Aenean lacinia bibendum nulla sed consectetur. Cum sociis natoque penatibus et magnis dis parturient montes, nascetur ridiculus mus. Donec ullamcorper nulla non metus auctor fringilla. Maecenas faucibus mollis interdum. Etiam porta sem malesuada magna mollis euismod.\",\"usage\":{\"executed\":\"35\",\"successful\":\"27\",\"positive\":\"12\",\"negative\":\"3\"}},{\"uuid\":\"uuid456\",\"text\":\"pBca 1256\",\"author\":{\"uuid\":\"uuid_author_456\",\"name\":\"Chris Anderson\",\"email\":\"chris@anderson.com\",\"biography\":\"This is a biography about Chris Anderson. It's different than Joe's... It's a little longer. Yada yada yada. Here's some latin. It should get truncated on the server or we could write our own directive to handle truncating (easy). Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.\"},\"description\":\"Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.\",\"usage\":{\"executed\":\"8\",\"successful\":\"8\",\"positive\":\"6\",\"negative\":\"0\"}}]}}");
            return msg.getJSONObject("data").toString();
        } catch (JSONException ex) {
            ex.printStackTrace();
            return null;
        }
        
    }
    
    //JCA:  as 0f 6/6/2013 submit seems to work
    public final void submit(String userText) {
        say(userText, "muted", true);
        if (!mind.runCommand(userText)) {
//            disambiguate(userText);  //JCA:  temporarily disabled for testing, also not fully hooked up
            say("Clotho was unable to satisfy that request", "text-error");
        } else {
            completer.put(userText);
        }
    }
    
    public final void learn(String nativeCmd, String jsCmd) {
        Interpreter.get().learnNative(nativeCmd, jsCmd);
    }
// </editor-fold> 

// <editor-fold defaultstate="collapsed" desc="Logging and Messaging"> 
    //JCA:  as 0f 6/9/2013 say seems to work
    public final void say(String message) {
        say(message, "text");
    }
    
    /**
     * 
     * @param message
     * @param severity  "text-error", "text", "text-warning", "text-success"
     *  see search-directives.js for 'from', is client or server
     */
    //JCA:  as 0f 6/9/2013 say seems to work
    public final void say(String message, String severity) {
//        System.out.println("say has : " + message);
        say(message, severity, false);

    }
    
    private final void say(String message, String severity, boolean isUser) {
        try {
            String source = "server";
            if(isUser) {
                source = "client";
            }
            JSONObject msg = new JSONObject();
                JSONObject data = new JSONObject();
                data.put("text", message);
                data.put("from", source);
                data.put("class", severity);
                data.put("timestamp", new Date().getTime());
            msg.put("channel", "say");
            msg.put("data", data);
            Router.get().sendMessage(mind.getClientConnection(), msg);
        } catch(Exception err) {
            err.printStackTrace();
        }
    }
    
    //JCA:  Java side looks fine, but client code crashes browser
    //clotho.alert("this is an alert!");
    public final void alert(String message) {
        System.out.println("alert has : " + message);

        try {
            JSONObject msg = new JSONObject();
            msg.put("channel", "alert");
            msg.put("data", message);
            Router.get().sendMessage(mind.getClientConnection(), msg);
        } catch(Exception err) {
            err.printStackTrace();
        }
    }
    
    
    //JCA:  This runs, and the message goes to the console.log spot.
    //clotho.log("I did some minipreps today");
    public final void log(String message) {
        System.out.println("log has : " + message);

        try {
            JSONObject msg = new JSONObject();
            msg.put("channel", "log");
            msg.put("data", message);
            Router.get().sendMessage(mind.getClientConnection(), msg);
        } catch(Exception err) {
            err.printStackTrace();
        }
    }
    
    //Make note of this message in my notebook
    public final void note(String message) {
        System.out.println("I need to put this in your notebook, but i'm not implemented");
        //These have the same structure as a say, but they are stored.
        
        say("I've stored your note (but not really): " + message);
    }
    
// </editor-fold> 

// <editor-fold defaultstate="collapsed" desc="Data Manipulation"> 

    //clotho.get("51b362f15076024ba019a642");  for a trail
    //clotho.get("51b29e7450765ced18af0d33");
    //JCA:  the object is requested, println'd, and collected in the clientside collector 6/8/2013
    //JCA:  result is also injected into the Mind's scriptengine as a proper object
    public final JSONObject get(String uuid) {

        try {
            //Pull the object from the db
            ObjBase datum = Collector.get().getObjBase(uuid.trim());
            Sharable obj = (Sharable) datum;
            JSONObject result = obj.toJSON();
            say("I found " + obj.getName(), "text-success");
            say("It has the value" + result.toString(), "text-success");
            
            //If the requestor is a client (not implemented) push the data to the client
            if(true) {
                JSONObject msg = makeCollect(obj);
                Router.get().sendMessage(mind.getClientConnection(), msg);
            }
            
            return result;
        } catch (Exception e) {
            say("Error retrieving " + uuid, "text-error");
        }
        
        return new JSONObject();
    }

    public final String set(String value) {
        try {
            //Pull the object from the db
            JSONObject newval = new JSONObject(value);
            
            //Check that a UUID has been provided
            String uuid = null;
            try {
                uuid = newval.getString("id");
            } catch(Exception err) {
                say("The arguments lack an 'id' field. Clotho does not know what object to alter");
                return null;
            }
            
            //Grab the object to be altered
            Sharable obj = null;
            try {
                ObjBase datum = Collector.get().getObjBase(uuid);
                obj = (Sharable) datum;
            } catch(Exception err) {
                say("No object with this id exists for Clotho to modify");
                return null;
            }
                    
                
            //Check that the current user has write access
            if(!Authenticator.get().hasWriteAccess(getPerson(), obj)) {
                say("The current user does not have write access for " + obj.getUUID().toString());
                return null;
            }
            
            //Complete the new data
            JSONObject existing = obj.toJSON();
            Iterator iterator = newval.keys();
            while(iterator.hasNext()) {
                String key = (String) iterator.next();
                Object val = newval.get(key);
                existing.put(key, val);
            }
            
            //Validate the data
            System.out.println("Stephanie needs to add set and create validation");
            if(false) {
                say("The data you wish to create did not pass validation.  No object was created" , "text-error");
                return null;
            }

            //Try to create the object and clobber the uuid
            String resultId = Persistor.get().save(existing);
            if(!resultId.equals(uuid)) {
                System.out.println("!resultId.equals(uuid) This should never happen");
                throw new Exception();
            }

            //Add the object to the Collector
            ObjBase object = Collector.get().temporaryRefetchMethod(uuid);

            //Contact the user to notify them that they modified an object
            say("You successfully modified: " + object.getName() + " with UUID: " + object.getUUID().toString(), "text-success");

            //Relay the data change to listening clients
            System.out.println("Ernst, this needs to be implemented.  Push object via pubsub.");

            //Return the modified data to the calling script
            return object.toString();
        } catch (Exception e) {
            say("Error setting " + value.toString(), "text-error");
            return null;
        }
    }

    public final String create(String json) {

    	try {
            //Determine whether the currently logged in person has permission to create
            if(!Authenticator.get().hasCreateAccess(getPerson())) {
                say("The current user does not have write access for this domain");
                return null;
            }

            //Construct the JSONObject
            JSONObject newval = new JSONObject(json);
            
            //Confirm that there is no pre-existing object with this uuid
            if(newval.has("id")) {
                String uuid = newval.getString("id");
                ObjBase datum = Collector.get().getObjBase(uuid);
                if(datum!=null) {
                    say("An object with the uuid " + uuid + " already exists.  No object was created." , "text-error");
                    return null;
                }
            }
            
            //Validate the data
            System.out.println("Stephanie needs to add set and create validation");
            if(false) {
                say("The data you wish to create did not pass validation.  No object was created" , "text-error");
                return null;
            }
            
            //Try to create the object and clobber the uuid
            JSONObject obj = new JSONObject(json);
            String uuidRes = Persistor.get().save(obj);
            
            //Add the object to the Collector and return its json
            ObjBase object = Collector.get().getObjBase(uuidRes);

            //Relay the data change to listening clients
            System.out.println("Ernst, this needs to be implemented here too.  Push object via pubsub.");

            //Return the JSON of the new object as a String
            say("You successfully created: " + object.getName() + " with UUID: " + object.getUUID().toString(), "text-success");
            return object.toString();
    	} catch(Exception e) {
            e.printStackTrace();
            return null;
    	}
    }
    
    //JCA:  as of 6/9/2013 works
    public final void destroy(String sharableId) {
        try {
            //Pull the object from the db and delete
            ObjBase datum = Collector.get().getObjBase(sharableId.trim());
            Sharable obj = (Sharable) datum;
            String name = obj.getName();
            String id = obj.getUUID().toString();
            Persistor.get().delete(obj);
            say("Sharable " + name + " with UUID " + id +  " has been destroyed", "text-success");
        } catch (Exception e) {
            say("Error destroying " + sharableId, "text-error");
        }
    }
    
    //clotho.query({"city","Townsville"});
    public final String query(String args) {
        try {
            //Execute the query
            System.out.println("Stephanie:  extract the query out more cleanly and define the scope of args (ssAPI.query())");
            JSONObject json = new JSONObject(args);
            HashMap<String, String> query = new HashMap<String, String>();
            Iterator iterator = json.keys();
            while(iterator.hasNext()) {
                String key = (String) iterator.next();
                query.put(key, (String) json.get(key));
            }

            //Relay the query to Persistor and return the hits
            System.out.println("Stephanie:  the persistor method should take an int for limiting the number of hits");
            //List<ObjBase> objs = Persistor.get().get(query, 100); //Should be like this
            List<ObjBase> objs = Persistor.get().get(query);
            
            //Wrap up results in a JSONArray
            JSONArray out = new JSONArray();
            for(ObjBase obj : objs) {
                try {
                    Sharable shar = (Sharable) obj;
                    //TODO: if current user has access to shar, otherwise ignore it
                    if(!Authenticator.get().hasReadAccess(getPerson(), shar)) {
                        continue;
                    }
                    JSONObject json2 = new JSONObject(shar.toString());
                    if(json2==null) {
                        continue;
                    }
                    out.put(json2);
                } catch(Exception err) {
                    err.printStackTrace();
                }
            }

            //Push results to the client (needs to move to context-specific)
            System.out.println("JCA:  you need to move this out to Router");
            for(ObjBase obj : objs) {
                try {
                    Sharable shar = (Sharable) obj;
                    if(shar==null) {
                        continue;
                    }
                    JSONObject msg = makeCollect(shar);
                    Router.get().sendMessage(mind.getClientConnection(), msg);
                } catch(Exception err) {
                    err.printStackTrace();
                }
            }
            
            //Say something to client in response
            System.out.println("JCA:  this too belongs only in Router");
            StringBuilder sb = new StringBuilder();
            for(ObjBase obj : objs) {
                sb.append("\n");
                sb.append(obj.getUUID().toString());
            }
            say("Clotho found " + out.length()+ " Sharables that satisfy your query: " + sb.toString(),  "text-success");
            return out.toString();
        } catch (Exception e) {
            say("Error querying ",  "text-error");
            return "[]";
        }
    }
// </editor-fold> 
    
// <editor-fold defaultstate="collapsed" desc="Execution"> 

    /**
     * Invoke the hard-coded editor appropriate for this particular object
     * *For an Instance, invoke the view set as the default view for this
     * type of data.  That's one option.  The other is to build up an editor
     * programmatically from the components.  Not sure.  Actually, I like this
     * second option best.  Maybe.  I don't know.  Second option is more secure.
     * 
     * Both options need to be available, I think.
     * 
     * @param sharableRef 
     */
    public final void edit(String sharableRef) {
    	/**
        try {
            Sharable sharable = (Sharable) resolveToObjBase(sharableRef);

            switch (sharable.type()) {
                case SCHEMA:
                    Schema schema = (Schema) sharable;
                    //Pop up the page, push in data
                    return;
                case INSTANCE:
                    return;
                default:
                    return;
            }
            //Pop up a new window of 
        } catch (Exception err) {
            err.printStackTrace();
        }
        **/
    }

    public final void run(String asstRef, String jsonArgs) {
    	/***
        Logger.log(Logger.Level.INFO, "received: " + asstRef + "   " + jsonArgs);
        try {
            Function assistant = (Function) resolveToObjBase(asstRef);
            JSONObject json = resolveToExecutionArguments(assistant,
                    resolveToJSONObject(jsonArgs));
            Logger.log(Logger.Level.INFO,
                    "Got Asst " + assistant.getName() + " w/ args " + json.toString());
            Executor.get().run(null, assistant, json, callback);
        } catch (JSONException e) {
            Logger.log(Logger.Level.WARN, "", e);
        }
        ***/    	
    }
    
    
    /**
     * Relay method for receiving a show command
     * @param viewRef
     * @param sharables
     * @param sloppy_args 
     */
    public final void show(String viewId, String sharables, String position) {
    	/***
        try {

            //IF THE PROCESS WAS TIED TO A COMMAND-INITIATED PROCESS THERE WOULD BE A PARENT DOO THAT NEEDS TO BE FOUND.
            Doo parentDoo = Hopper.get().extract(null);
            ShowDoo doo = new ShowDoo(parentDoo);

            //Gather up referenced objects
            View view = (View) resolveToObjBase(viewId);
            List<Sharable> shareList = resolveToSharables(sharables);

            //Create the widget and put it on its page
            Page targetPage = mind.getConfig().getPage(socket_id);
            Widget widget = new Widget(targetPage, view);
            targetPage.addWidget(widget);

            //Just for record keeping
            doo.viewId = view.getId();
            doo.widgetId = widget.getId();
            doo.collectSharables(shareList);

            //Create all the commands for the client
            doo.commandMessageArray = new JSONArray();
            doo.commandMessageArray.put(makeCollect(shareList));
            doo.commandMessageArray.put(makeShowWidget(widget, position));
            doo.commandMessageArray.put(makeUpdate(widget, shareList, socket_id));
            doo.commandMessageArray.put(makeCallback(doo.getId(), socket_id));

            //Put the doo into the hopper to await a callback
            Hopper.get().add(doo);

            //Send the commands
            Communicator.get().sendClientMessage(socket_id,
                    SendChannels.commandList,
                    doo.commandMessageArray.toString());

            //Save everything whose state was changed  //JCA:  THIS NEEDS A CALLBACK/FAILURE RESPONSE THAT REVERTS THIS (EVENTUALLY)
            Persistor.get().persistObjBase(mind);
        } catch (Exception e) {
            Logger.log(Logger.Level.WARN, "", e);
            e.printStackTrace();
            //REVERT THE MIND
            //SEND CLIENT MESSAGE TELLING THAT IT FAILED
        }
        ***/
    }

    public final void startTrail(String trailRef) {
    	/***
        try {
            //Grab a Doo if something is awaiting one
            Doo parentDoo = Hopper.get().extract(null);
            Trail trail = (Trail) resolveToObjBase(trailRef);
            Person student = this.getPerson();

            Proctor.get().initiateTrail(student, trail, parentDoo);

        } catch (Exception e) {
            Logger.log(Logger.Level.WARN, "", e);
            e.printStackTrace();
        }
        ***/
    }
    
// </editor-fold> 

    // <editor-fold defaultstate="collapsed" desc="Setup"> 
 
    public ServerSideAPI(ClothoWebSocket socket) {
    	this.socket = socket;
    	this.mind = null;
    }
    public ServerSideAPI(Mind mind) {
        this.mind = mind;
    }
    
    public void setSocket(ClothoWebSocket socket) {
    	this.socket = socket;
    }
    public void setMind(Mind mind) {
    	this.mind = mind;
    }
// </editor-fold>

//    public final void test() {
//        try {
//            System.out.println("Test has been invoked");
//            JSONObject trail = new JSONObject("{\"uuid\":\"trail_123\",\"title\":\"Biosafety Module\",\"author\":\"UC Berkeley\",\"description\":\"<blockquote><p>This is a module on Biosafety. You'll learn about Cras sit amet nibh libero, in gravida nulla. Nulla vel metus scelerisque ante sollicitudin commodo. Cras purus odio, vestibulum in vulputate at, tempus viverra turpis.</p><p><small>Maxwell Bates</small></p></blockquote>\",\"contents\":[{\"module_title\":\"The Basics\",\"pavers\":[{\"paver_title\":\"Introduction\",\"type\":\"template\",\"template\":\"/app/partials/trail_123_1.html\"}]},{\"module_title\":\"Biosafety Levels\",\"pavers\":[{\"paver_title\":\"Introduction\"},{\"paver_title\":\"Biosafety Level 1-2\"},{\"paver_title\":\"Biosafety Level 3-4\"}]},{\"module_title\":\"Assessment\",\"pavers\":[{\"paver_title\":\"Review\"}],\"assessment\":[{\"type\":\"quiz\",\"title\":\"Final Quiz\"}]}]}");
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
     * @return 
     */
    private Person getPerson() {
    	return mind.getPerson();
    	/**
        if (person != null) {
            return person;
        }
        return Person.getById(mind.getPersonId());
        **/    	
    }

    public final void newPage() throws Exception {
    	/***
        try {
            JSONArray commandMessageArray = new JSONArray();
            commandMessageArray.put(makeNewPage(UUID.randomUUID().toString()));

//            
//            
//            //Put the doo into the hopper to await a callback
//            Hopper.get().add(doo);

            //Send the commands
            Communicator.get().sendClientMessage(socket_id,
                    SendChannels.commandList,
                    commandMessageArray.toString());

            //Save everything whose state was changed  //JCA:  THIS NEEDS A CALLBACK/FAILURE RESPONSE THAT REVERTS THIS (EVENTUALLY)
            Persistor.get().persistObjBase(mind);
        } catch (Exception e) {
            Logger.log(Logger.Level.WARN, "", e);
            e.printStackTrace();
            //REVERT THE MIND
            //SEND CLIENT MESSAGE TELLING THAT IT FAILED
        }
		***/
    }


    /**
    //if a command has a callback, handle it here
    //minimally removes the doo from the dooHopper   
    //jca:  i don't know that i like this being in the api like this.  it really
    //should be on a more secure channel because this all deals with Doos, which
    //need to be accurately tracked.  if you can directly call this callback from
    //the api, you can just type it into the commmand bar and lie to clotho, which
    //is no good.  so, yeah, it shouldn't be implemented like this, but we can fix that
    //in the future when we pick through the client-->server communications standard.
    //also, this is a mind-agnostic thing -- this is what gets called just to confirm
    //that a process has completed, and everything should probably call this.  The doo
    //hopper maybe should even be an Aspect with this being relayed there.  Yeah.  That's
    //it.*/
    public final void notify(String dooID) {
    	/**
        System.out.println("calling back!");
        //Grab the Doo and perhaps perform any commands embedded in there

        //Terminate it from the Hopper
        Hopper.get().terminate(dooID);
        **/
    }
    
    // <editor-fold defaultstate="collapsed" desc="Command Factories"> 


    /**
     * Instruct the client to receive a List of Sharables and log
     * those objects into the clientside collector
     * 
     * @param sharables
     * @return
     * @throws Exception 
     */
    public static JSONObject makeCollect(Sharable sharable) throws Exception {
        JSONObject msg = new JSONObject();
            JSONObject data = new JSONObject();
            data.put("uuid", sharable.getUUID());
            data.put("type", "json");
            data.put("model", sharable.toJSON());
            data.put("isURL", "false");
            
        msg.put("data", data);
        msg.put("channel", "collect");

//        System.out.println("Make Collect has:\n " + msg.toString(2));
        return msg;
    }

    /**
     * Instantiate a new workspace page.  Only the Proctor can construct
     * arguments for TRAILS or specific editors.  Security thing.
     * @return
     * @throws Exception 
     */
    public static JSONObject makeNewPage(String page_id) throws Exception {
        JSONObject out = new JSONObject();
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
    public JSONObject makeShowWidget(Widget widget, String position) throws Exception {
        JSONObject out = new JSONObject();
        JSONObject positioninfo = resolveToJSONObject(position);
        View view = widget.getView();

        //Load the html and js scripts into a JSONObject, insert the uuid for the widget
        String widgetId = widget.getId();
        System.out.println("graphics script: "+view.getGraphicsScript());
        String html = replaceWidgetId(view.getGraphicsScript(), widgetId);
        String onshow = replaceWidgetId(view.getOnShowScript(), widgetId);

        out.put("widget_id", widgetId);
        out.put("content", html);
        out.put("on_show", onshow);
        out.put("parent_widget_id", "widget_space_root");
        out.put("command", "showWidget");

        return out;
    }

    /**
     * Update instructs the client to call update on a widget
     * The data that goes into the widget, however, is sent in a Collect command
     * 
     * @param view
     * @return
     * @throws Exception 
     */
    public static JSONObject makeUpdate(Widget widget, List<Sharable> sharables, String socket_id) throws Exception {
        JSONObject updateCmd = new JSONObject();
        //Createhte token:uuid map that associates the piece of data with its token used by the widget
        List<ClothoField> fields = widget.getView().getInputArguments();
        JSONObject inputs = new JSONObject();
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

    public static JSONObject makeCallback(String dooId, String socket_id) throws Exception {
        JSONObject callBackCmd = new JSONObject();
        callBackCmd.put("socketID", socket_id);
        callBackCmd.put("dooID", dooId);
        callBackCmd.put("command", "callback");

        return callBackCmd;
    }



    /***                                ***\
     ********  Utility methods  ********
    \*                                  */
    //XXX: needs updating
    private JSONObject resolveToExecutionArguments(Function asst, JSONObject input) throws JSONException {
        //Do type validation of the arguments, not call candooit
        //List<ClothoField> inputTypes = asst.;
        //return argumentsRelay(inputTypes, input);
        throw new UnsupportedOperationException();
    }

    private static JSONObject argumentsRelay(List<ClothoField> inputTypes,
            JSONObject input) throws JSONException {
        //For each token of this schema, see if the input object has the token
        /* TODO: is this right? */
        for (ClothoField field : inputTypes) {
            String token = field.getName();
            if (!input.has(token)) {
                break;
            }
            return input;
        }

        //If it got this far, there is something wrong with the input arguments (not wrapped properly)

        //See if it is just a single-token thing with a single-object argument needs to be wrapped up
        try {
            JSONObject out = new JSONObject();
            for (ClothoField field : inputTypes) {
                String token = field.getName();
                out.put(token, input);
            }
            return out;
        } catch (Exception err2) {
        }
        return null;
    }

    /**
     * Inputs a sloppy input of a datum reference and finds the 'correct' server-side
     * object they are talking about
     * @param schemaRef
     * @return 
     */
    private ObjBase resolveToObjBase(Object datumRef) {
        ObjBase out = null;
        Logger.log(Logger.Level.INFO, "resolveToObjBase has " + datumRef);

        //If it's a ObjBase object
        try {
            ObjBase aschmema = (ObjBase) datumRef;
            ObjectId uuid = aschmema.getUUID();

            //Fetch the real version of the ObjBase
            out = Collector.get().getObjBase(uuid.toString());
            if (out != null) {
                return out;
            }
        } catch (Exception e) {
            Logger.log(Logger.Level.INFO, "It wasn't a ObjBase object");
        }

        //If it's a String
        try {
            String str = (String) datumRef;

            //See if it's a uuid String
            out = Collector.get().getObjBase(str);
            if (out != null) {
                return out;
            }

            //See if it's a json string representation of a ObjBase
            JSONTokener tokener = new JSONTokener(str);
            JSONObject jsonobj = new JSONObject(tokener);

            String uuid = jsonobj.getString("uuid");
            out = Collector.get().getObjBase(uuid);
            if (out != null) {
                return out;
            }

            //See if it's the proper value of a NameField of some ObjBase
            //THIS INVOLVES COLLATOR FUNCTIONALITY

        } catch (Exception e) {
            Logger.log(Logger.Level.INFO, "It wasn't a string of json");
        }

        return out;
    }

    /**
     * Inputs a sloppy argument of some JSON -- like whether the user
     * put in a String of it, Java object representation, or whatever
     * @param jsonArgs
     * @return 
     */
    public final JSONObject resolveToJSONObject(Object jsonArgs) {
        System.out.println(jsonArgs.getClass().toString());
        //If it's a JSONObject
        try {
            JSONObject out = (JSONObject) jsonArgs;
            return out;
        } catch (Exception e) {
//            Logger.log(Logger.Level.INFO, "It wasn't a jsonobject");
        }
        
        //If it's a String representation of json
        try {
            String str = (String) jsonArgs;
            JSONTokener tokener = new JSONTokener(str);
            JSONObject out = new JSONObject(tokener);
            return out;
        } catch (Exception err) {
//            Logger.log(Logger.Level.INFO, "It wasn't a string of json");
        }

        //If it's a String uuid
        try {
            String uuid = (String) jsonArgs;
            ObjBase d = Collector.get().getObjBase(uuid);
            Sharable share = (Sharable) d;
            return share.toJSON();
        } catch (Exception err) {
//            Logger.log(Logger.Level.INFO, "It wasn't a string of uuid");
        }

        
        //If it's some other kind of object
        Logger.log(Logger.Level.WARN, String.valueOf(jsonArgs));
        return null;
        
    }

    /**
     * Resolves a sloppy list of sharables from a String into a List of
     * occupied uuids.  Null or non-existent uuids return null
     * 
     * @param sharables
     * @return
     * @throws Exception 
     */
    public static final List<Sharable> resolveToSharables(String sharables) throws Exception {
        ArrayList<Sharable> out = new ArrayList<Sharable>();
        JSONArray array = resolveToJSONArray(sharables);
        for (int i = 0; i < array.length(); i++) {
            String uuid = array.getString(i);
            Sharable sharable = (Sharable) Collector.get().getObjBase(uuid);
            out.add(sharable);
        }
        return out;
    }

    /**
     * Inputs a sloppy argument of some JSON -- like whether the user
     * put in a String of it, Java object representation, or whatever
     * @param jsonArgs
     * @return 
     */
    public static final JSONArray resolveToJSONArray(Object jsonArgs) {
        //If it's a JSONArray
        try {
            JSONArray out = (JSONArray) jsonArgs;
            return out;
        } catch (Exception err) {
        }

        //If it's a String representation of a JSONArray
        try {
            String str = (String) jsonArgs;
            JSONTokener tokener = new JSONTokener(str);
            JSONArray out = new JSONArray(tokener);
            return out;
        } catch (Exception err) {
        }

        //If it's a single String like one uuid
        try {
            String str = (String) jsonArgs;
            JSONArray out = new JSONArray();
            out.put(str);
            return out;
        } catch (Exception err) {
        }

        //If it's some other kind of object
        Logger.log(Logger.Level.WARN, "The jsonArgs could not be resolved");
        return null;
    }

    /**
     * Replace the _widget_id phrases with the actual uuid
     * @param script
     * @param widgetIdPrefix
     * @return 
     */
    public static final String replaceWidgetId(String script, String widgetIdPrefix) {
        try {
            return XMLParser.addPrefixToTagAttribute(script, "id", widgetIdPrefix);
        } catch (Exception ex) {
        	Logger.log(Logger.Level.FATAL, ex);
        }
        return null;
    }

        private void disambiguate(String nativeCmd) {
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

                if (mind.runCommand(assoc_cmd)) {
                    // learn to associate nativeCmd with assoc_cmd
                } else {
                    // tell the user that native cmd failed, so the overall action requested was aborted.
                    //Commicator.get().say("sorry man...no dice");
                }
                //IT SHOULD RETURN THE SEARCH RESULTS TO THE CLIENT, THEN LET THEM EXECUTE, BUT SINCE WE DON'T HAVE THAT NOW, IT WILL JUST RUN IT
                //I SEE A POTENTIAL ISSUE HERE OF NOT HAVING ACCESS TO THE MIND IN QUESTION
            }
        }
    
    // </editor-fold> 
    
        
    /********  VARIABLES  ********/
    //private transient Person person;
    private static final AutoComplete completer = new AutoComplete();
    private String socket_id;
    private ClothoWebSocket socket;
    private Mind mind;
}

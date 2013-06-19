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

import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import javax.swing.JOptionPane;
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
import org.clothocad.core.util.FileUtils;
import org.clothocad.model.Person;
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

	private Persistor persistor;
	
	public ServerSideAPI() {
		try {
			this.persistor = new Persistor(new MongoDBConnection());
		} catch(Exception e) {
			e.printStackTrace();
		}
		
		// here, we could also create a pool of persistors
	}

// <editor-fold defaultstate="collapsed" desc="Human Interaction">      
    //JCA:  as of 6/6/2013 autcomplete works.  Wordlist is not persisted, but the completer does learn submitted phrases.
    public final String autocomplete(String rawMsg) {
        String userText = null;
        try {
            JSONObject json = new JSONObject(rawMsg);
            userText = json.getString("query");
        } catch(Exception err) {
            userText = rawMsg;
        }
        
        try {
            ArrayList<String> completions = completer.getCompletions(userText);
            JSONArray data = new JSONArray(); //The list of JSONObjects each for an autocomplete
            for(String str : completions) {
                JSONObject obj = new JSONObject();
                obj.put("text", str);
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
    public final void submit(String command) {
        //Resolve the arguments to a command string
        String message = null;
        try {
            JSONObject obj = new JSONObject(command);
            
            //Extract the messsage
            if(obj.has("query")) {
                message = obj.getString("query");
            } else if(obj.has("message")) {
                message = obj.getString("message");
            } else if(obj.has("value")) {
                message = obj.getString("value");
             } else if(obj.has("msg")) {
                message = obj.getString("msg");
            }
        } catch(Exception notAJSON) {
            message = command;
        }
        
        
        say(message, "muted", null, true);
        
        if (!mind.runCommand(message)) {
            disambiguate(message);  //JCA:  temporarily disabled for testing, also not fully hooked up
            return;
        }
        
        //If the command successfully executed, it gets retained
        mind.addLastCommand(message);
    }
    
    //clotho.learn("test1", "clotho.say('Hi there');");
    
    public final void learn(String nativeCmd, String jsCmd) {
        JSONObject json = null;
        
        //See if the jsCmd is actually an execution statement
        try {
            json = new JSONObject(jsCmd);
        } catch(Exception err) {
            
            //See if it's empty, meaning use the last command they issued
            if(jsCmd==null) {
                json = mind.getLastCommands().get(0);
                if(json==null) {
                    System.out.println("This is something that will happen if you haven't issued a previous command, it has nothing to learn");
                    return;
                }
                
            //Otherwise it is a js stqatement
            } else {
                try {
                    json = new JSONObject();
                    json.put("text", jsCmd);
                    json.put("type", "phrase");
                } catch (JSONException ex) {
                    System.out.println("This should never happen");
                    ex.printStackTrace();
                    return;
                }
            }
        }
        
        Interpreter.get().learnNative(nativeCmd, json);
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
        say("The mind has been cleared", "text-success");
    }
// <editor-fold defaultstate="collapsed" desc="Logging and Messaging"> 
    //JCA:  as 0f 6/9/2013 say seems to work
    public final void say(String message) {
        say(message, "text", null, false);
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
        say(message, severity, null, false);

    }
    
    public final void say(String message, String severity, String recipients) {
//        System.out.println("say has : " + message);
        say(message, severity, recipients, false);

    }
    
    private final void say(String rawMsg, String severity, String recipients, boolean isUser) {
        String message = null;
        //If rawMsg is an object, parse out the arguments
        try {
            JSONObject obj = new JSONObject(rawMsg);
            
            //Extract the messsage
            if(obj.has("message")) {
                message = obj.getString("message");
            } else if(obj.has("value")) {
                message = obj.getString("value");
            } else if(obj.has("msg")) {
                message = obj.getString("msg");
            }
            
            //If the message is still null, println the whole object
            if(message==null) {
                message = obj.toString();
                
            //Otherwise it extracted out a message, so keep parsing
            } else {

                //Extract the severity
                if(obj.has("severity")) {
                    severity = obj.getString("severity");
                }

                if(obj.has("recipients")) {
                    recipients = obj.getString("recipients");
                }
            }
        } catch(Exception notAJSON) {
            message = rawMsg;
        }
        
        //Resolve the recipients
        List<Sharable> listUsers = null;
        if(recipients !=null) {
            listUsers = resolveToExistentSharablesList(recipients);
        }
        
        //Assign the say as server or client speaking (not user defined, is logic-based);
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
    public final String get(String sharableRef) {
        String out = null;
        try {
            //Resolve the arguments and retrieve
            Sharable existing = resolveToExistentSharable(sharableRef);
            if(existing==null) {
                throw new Exception();
            }
            //Check that the user has permission
            if(!Authenticator.get().hasReadAccess(getPerson(), existing)) {
                say("The current user does not have read access for " + sharableRef, "text-error");
                System.out.println("JCA:  this error msg should be the same as the previous.  The object should appear to not exist if denied");
                return null;
            }
            
            //If the requestor is a client (not implemented) push the data to the client
            if(true) {
                JSONObject msg = makeCollect(existing);
                Router.get().sendMessage(mind.getClientConnection(), msg);
                
                //Associate the client with this sharable in pubsub
                System.out.println("Ernst:  this is the registering of pubsub");
                Router.get().register(mind.getClientConnection(), existing);
            }
            
            //Notify the user that they retrieved the data
            System.out.println("JCA:  this should be moved specifically to search bar responses");
            
            out = existing.toString();
            say("I retrieved the Sharable " + out, "text-success");

        } catch (Exception e) {
            //Start of fudgy retieval from filesystem
            try {
                String uuid = this.relayResolveStringToUUID(sharableRef);
                System.out.println("Stephanie:  JCA hack to deal with filesystem-persisted models.  Should move to db.");
                out = FileUtils.readFile("clotho3-web/models/" + uuid + ".json");
                JSONObject obj = new JSONObject(out);
                
                JSONObject msg = new JSONObject();
                    JSONObject data = new JSONObject();
                    data.put("uuid", sharableRef);
                    data.put("type", "json");
                    data.put("model", obj);
                    data.put("isURL", "false");

                msg.put("data", data);
                msg.put("channel", "collect");
                Router.get().sendMessage(mind.getClientConnection(), msg);
                return out;
                
            } catch(Exception err) {
                say("Error getting from filesystem " + sharableRef, "text-error");
            }
            //End of fudgy retieval from filesystem
                            
                
                
            //Start of super-fudgy short-circuit
            
                try {
                    if(sharableRef.equals("org.clothocad.model.Institution")) {
                        System.out.println("Stephanie, I need for the Collector.getObjBase request to return the Schema, and then retrieve the JSONObject representation");
                        JSONObject jsonSchema = new JSONObject("{\"schema\":[{\"name\":\"name\",\"readable\":\"Display Name\",\"type\":\"text\",\"placeholder\":\"Your Name\",\"required\":true},{\"name\":\"id\",\"readable\":\"Id\",\"type\":\"text\",\"placeholder\":\"Your Name\",\"required\":true},{\"name\":\"city\",\"readable\":\"City\",\"type\":\"text\",\"placeholder\":\"Your City\",\"required\":true},{\"name\":\"state\",\"readable\":\"State\",\"type\":\"text\",\"placeholder\":\"Your State\",\"required\":true},{\"name\":\"country\",\"readable\":\"Country\",\"type\":\"text\",\"placeholder\":\"Your Country\",\"required\":true},],\"custom\":{}}");
                        JSONObject msg = new JSONObject();
                            JSONObject data = new JSONObject();
                            data.put("uuid", "org.clothocad.model.Institution");
                            data.put("type", "json");
                            data.put("model", jsonSchema);
                            data.put("isURL", "false");

                        msg.put("data", data);
                        msg.put("channel", "collect");
                        Router.get().sendMessage(mind.getClientConnection(), msg);
                        return jsonSchema.toString();
                    }
                } catch (JSONException ex) {
                    java.util.logging.Logger.getLogger(ServerSideAPI.class.getName()).log(Level.SEVERE, null, ex);
                }
            
            //End of super-fudgy short-circuit

            
            
            
            
            

        }
        
        return out;
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

            }
            
            //To deal with nested data (ideally this wouldn't be here, but it is)
            if(uuid==null) {
                try {
                    newval = newval.getJSONObject("data");
                    uuid = newval.getString("id");
                } catch(Exception err) {
                say("The arguments lack an 'id' field. Clotho does not know what object to alter");
                return null;
                }
            }
            //End: deal with nested data
            
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
            
            //Confirm that the new data is different than the old data
            JSONObject original = obj.toJSON();
            if(original.toString().equals(existing.toString())) {
                say("The data was unmodified." , "text-warning");
                return original.toString();
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
            Router.get().publish((Sharable) object);

            //Return the modified data to the calling script
            return object.toString();
        } catch (Exception e) {
            say("Error setting " + value.toString(), "text-error");
            return null;
        }
    }

    public final void set(JSONObject data) {
    	
    }
    
    /**
    public final String create(JSONObject json) {
    
    	//System.out.println("[ServerSideAPI.create] " + json);
    	
    	try {
    		
    		String uuid = json.getString("uuid");

    		JSONObject model = json.getJSONObject("model");
    		

    		/*** dynamic class loading  
    		// this needs to be improved somehow...
    		String sClass = json.getString("className");
    		Class<?> c = Class.forName(sClass);
    		
    		// how can I figure out which class this is?
    		NucSeq nucseq = (NucSeq)c.newInstance();

    		// now, we need to set the sequence of this nucseq object
    		nucseq.initiateNucSeq(model.getString("sequence"));
    		    		    		    		    	
    		// can/should we create an ObjBase object here from 
    		// the provided JSON information??
    				
    		// here we need to forward the model to the Persistor
    		if(null != this.persistor) {
    			persistor.save(nucseq);
        		return nucseq.getUUID().toString();
    		} else {
    			System.err.println("CRAP!");
    		}
    		
    		return new String();
    		
    		
    	} catch(Exception e) {
    		e.printStackTrace();
    		return e.getMessage();
    	}
    }
    **/
    
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
            if(uuidRes==null) {
                say("The object could not be persisted during create ", "text-error");
                return null;
            }
            
            //Add the object to the Collector and return its json
            ObjBase object = Collector.get().getObjBase(uuidRes);
            if(object==null) {
                say("The object was created, but could not be retrieved ", "text-error");
                return null;
            }

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
            List<Sharable> existent = this.resolveToExistentSharablesList(sharableId);
            for(Sharable obj : existent) {
                //Check that the current user has write access
                if(!Authenticator.get().hasWriteAccess(getPerson(), obj)) {
                    say("The current user does not have write access for " + obj.getUUID().toString());
                    continue;
                }

                String name = obj.getName();
                String id = obj.getUUID().toString();
                Persistor.get().delete(obj);
                say("Sharable " + name + " with UUID " + id +  " has been destroyed", "text-success");
            }
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

                    if(!Authenticator.get().hasReadAccess(getPerson(), shar)) {
                        continue;
                    }
                    
                    JSONObject json2 = shar.toJSON();
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
    public final void show(String viewId, String sharables, String position, String recipients) {
        System.out.println("Show is called " + position);
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
            //Resolve the arguments and retrieve, this will also push refreshed data to client and register pubsub
            String existing = get(sharableRef);
            if(existing==null) {
                 say("Clotho was unable to resolve the arguments for edit", "text-error");
                 return;
            }
                
            try {
                JSONObject json = new JSONObject(existing);
                String uuid = json.getString("id");
                

                
                
                JSONObject msg = new JSONObject();
                    JSONObject data = new JSONObject("{\"template\":\"widget/dependencies/simple-template.html\",\"target\":\"body\",\"controller\":\"widget/dependencies/simple-controller.js\",\"dependencies\":[\"widget/dependencies/simple-service.js\"],\"styles\":{\"background-color\":\"#FF0000\"}}");
                    data.put("args", uuid);

                msg.put("data", data);
                msg.put("channel", "display_simple");
                
                
                Router.get().sendMessage(mind.getClientConnection(), msg);
            } catch (JSONException ex) {
                 say("Clotho was unable to invoke edit", "text-error");
                 return;
            }
            

    }
    
    public final void listen(String args) {
        say("not yet implemented", "text-error");
    }

    public final void unlisten(String data) {
        say("not yet implemented", "text-error");
    }

// </editor-fold> 

    // <editor-fold defaultstate="collapsed" desc="Setup"> 
 
    public ServerSideAPI(ClothoWebSocket socket) {
    	this.socket = socket;
    	this.mind = null;

    	try {
			this.persistor = new Persistor(new MongoDBConnection());
		} catch (Exception e) {
			e.printStackTrace();
		}
    }
    
    public ServerSideAPI(Mind mind) {
        this.mind = mind;
        
        try {
			this.persistor = new Persistor(new MongoDBConnection());
		} catch (Exception e) {
			e.printStackTrace();
		}
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
//        JSONObject out = new JSONObject();
//        JSONObject positioninfo = resolveToJSONObject(position);
//        View view = widget.getView();
//
//        //Load the html and js scripts into a JSONObject, insert the uuid for the widget
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
// </editor-fold> 
    
    // <editor-fold defaultstate="collapsed" desc="Utility Methods"> 

    private List<Sharable> resolveToExistentSharablesList(String sharableList) {
        ArrayList<String> uuidList = new ArrayList<String>();
        
        //If it's a JSON Array of stuff
        try {
            JSONArray ary = new JSONArray(sharableList);
            for(int i=0; i<ary.length(); i++) {
                String str = ary.getString(i);
                Sharable shar = resolveToExistentSharable(str);
                uuidList.add(shar.getUUID().toString());
            }
        } catch(Exception notJsonArray) {
            //If it's a single object
            try {
                JSONObject obj = new JSONObject(sharableList);
                String uuid = relayResolveJSONObjectToUUID(obj);
                uuidList.add(uuid);
            } catch(Exception notJsonObject) {
                try {
                    this.relayResolveStringToUUID(sharableList);
                } catch(Exception notStringRef) {
                    return null;
                }
            }
        }
        
        ArrayList<Sharable> out = new ArrayList<Sharable>();
        for(String uuid : uuidList) {
            ObjBase obj = Collector.get().getObjBase(uuid);
            Sharable shar = (Sharable) obj;
            out.add(shar);
        }
        return out;
    }
    
    private Sharable resolveToExistentSharable(String sharableRef) {
        //Extract UUIDs based on how they were expressed
        String uuid = null;
        
        try {
            JSONObject obj = new JSONObject(sharableRef);
            uuid = relayResolveJSONObjectToUUID(obj);
        } catch(Exception notJsonObject) {
            try {
                JSONArray ary = new JSONArray(sharableRef);
                uuid = relayResolveJSONArrayToFirstUUID(ary);
            } catch(Exception notJSONArray) {
                try {
                    uuid = relayResolveStringToUUID(sharableRef);
                } catch(Exception notStringArgs) {
                    return null;
                }
            }
        }
        
        //Pull the Sharables by ID, put in a List, and return them
        ObjBase obj = Collector.get().getObjBase(uuid);
        Sharable sharable = (Sharable) obj;
        return sharable;
    }
        
        private String relayResolveJSONArrayToFirstUUID(JSONArray ary) throws Exception {
            try {
                JSONObject json = ary.getJSONObject(0);
                return relayResolveJSONObjectToUUID(json);
            } catch(Exception notJSON) {}
            
            //Try to parse as a UUID in a list, if fails will throw exception
            String data = ary.getString(0);
            ObjectId oid = new ObjectId(data);
            return data;
        }
    
        private String relayResolveJSONObjectToUUID(JSONObject obj) throws Exception {
            //Scenario:  the object is a single JSON with an id
            if(obj.has("id")) {
                try {
                    String uuid = obj.getString("id");
                    return uuid;
                } catch(Exception err) {}
            }
            
            //Scenario:  they passed an object or array wrapped as 'data', 'value', 'args', 'obj', or 'object'
            try {
                return relayResolveToObjectWithToken(obj, "data");
            } catch(Exception err) {};
            try {
                return relayResolveToObjectWithToken(obj, "value");
            } catch(Exception err) {};
            try {
                return relayResolveToObjectWithToken(obj, "args");
            } catch(Exception err) {};
            try {
                return relayResolveToObjectWithToken(obj, "obj");
            } catch(Exception err) {};
            try {
                return relayResolveToObjectWithToken(obj, "object");
            } catch(Exception err) {};

            throw new Exception();
        }
        
            private String relayResolveToObjectWithToken(JSONObject obj, String token) throws Exception {
                if(obj.has(token)) {
                    try {
                        String data = obj.getString(token);
                        ObjectId oid = new ObjectId(data);
                        return data;
                    } catch(Exception err) {}
                    try {
                        JSONObject data = obj.getJSONObject(token);
                        String uuid = data.getString("id");
                        return uuid;
                    } catch(Exception err) {}
                    try {
                        JSONArray ary = obj.getJSONArray(token);
                        JSONObject data = ary.getJSONObject(0);
                        String uuid = data.getString("id");
                        return uuid;
                    } catch(Exception err) {}
                }
                throw new Exception();
            }
        
        //Scenario:  the argument did not parse into a JSONObject, it is either
        //the name of a Sharable, a namespaced token, or a direct UUID
        private String relayResolveStringToUUID(String str) throws Exception {
            //Try parsing as uuid
            try {
                ObjectId id = new ObjectId(str.trim());
                return id.toString();
            } catch(Exception err) {}
            
            //Check if it's in the Person's namespace
            String uuid = mind.pullUUIDFromNamespace(str);
            if(uuid!=null) {
                return uuid;
            }
            
            //Try a query by name and allow to fail
            try {
                HashMap map = new HashMap();
                map.put("name", str);
                List<ObjBase> result = Persistor.get().get(map);
                ObjBase first = result.get(0);
                return first.getUUID().toString();
            } catch(Exception err) {}
            
            return str;
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
        try {
            say("Entering disambiguation");
            Set<String> cmdResults = Interpreter.get().receiveNative(nativeCmd);

            if(cmdResults.isEmpty()) {
                say("No suggestions are available.", "text-warning");
                return;
            }
            
            JSONArray data = new JSONArray();
            for(String cmdResult : cmdResults) {
                try {
                    JSONObject oneSuggestion = new JSONObject(cmdResult);
                    data.put(oneSuggestion);
                } catch (JSONException ex) {
                    ex.printStackTrace();
                }
                say(cmdResult);
            }
            JSONObject msg = new JSONObject();
            msg.put("channel", "autocomplete");
            msg.put("data", data);
            Router.get().sendMessage(mind.getClientConnection(), msg);
        } catch (JSONException ex) {
            ex.printStackTrace();
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

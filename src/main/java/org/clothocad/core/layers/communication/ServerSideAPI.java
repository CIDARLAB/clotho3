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
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.UUID;
import javax.swing.JOptionPane;

import javax.xml.parsers.ParserConfigurationException;

import org.clothocad.core.aspects.Collector;
import org.clothocad.core.aspects.Interpreter.AutoComplete;
import org.clothocad.core.aspects.Logger;
import org.clothocad.core.aspects.Persistor;
import org.clothocad.core.aspects.Interpreter.Interpreter;
import org.clothocad.core.datums.Datum;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.Sharable;
import org.clothocad.core.datums.View;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.layers.communication.connection.ClientConnection;
import org.clothocad.core.layers.communication.connection.ws.ClothoWebSocket;
import org.clothocad.core.layers.communication.mind.Mind;
import org.clothocad.core.layers.communication.mind.Widget;
import org.clothocad.core.layers.persistence.mongodb.MongoDBConnection;
import org.clothocad.core.schema.Schema;
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
 * So, it's just a wrapper really.
 * 
 * @author John Christopher Anderson
 */
public final class ServerSideAPI {

    private String socket_id;
    private ClothoWebSocket socket;
    private Mind mind;
    
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
    
    /***                                ***\
     ********  Public API methods  ********
    \*                                  */
    

    public final void autocomplete(String userText) {
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
            msg.put("channel", "autcomplete");
            msg.put("data", data);
            Router.get().sendMessage(mind.getClientConnection(), msg);
        } catch(Exception err) {
            err.printStackTrace();
        }
    }
    
    //JCA:  as 0f 6/6/2013 submit seems to work
    public final void submit(String userText) {
        if (!mind.runCommand(userText)) {
            disambiguate(userText);
            completer.put(userText);
        }
    }

    //JCA:  as 0f 6/6/2013 say seems to work
    public final void say(String message) {
        System.out.println("say has : " + message);

        try {
            JSONObject msg = new JSONObject();
                JSONObject data = new JSONObject();
                data.put("text", message);
                data.put("from", "toBeWorkedOutLater");
                data.put("class", "text-error");
                data.put("timestamp", new Date().getTime());
            msg.put("channel", "say");
            msg.put("data", data);
            Router.get().sendMessage(mind.getClientConnection(), msg);
        } catch(Exception err) {
            err.printStackTrace();
        }
    }
    
    public final JSONObject get(String uuid) {
    	/**
        try {
            Datum datum = Collector.get().getDatum(uuid);
            Sharable obj = (Sharable) datum;
//            Logger.log(Logger.Level.INFO, "received sharable: " + uuid);
            JSONObject result = obj.toJSON();
            return result;
        } catch (ClassCastException e) {
            Logger.log(Logger.Level.WARN,
                    "could not find sharable: " + uuid,
                    e);
        }
        **/
        return new JSONObject();
    }

    public final void set(String sharableId, String newvalue) {
    	/**
        Sharable sharable = (Sharable) resolveToDatum(sharableId);
        JSONObject value = this.resolveToJSONObject(newvalue);
        if (sharable.set(value, this.getPerson(), new ShowDoo(null))) {
            //RETURN A COMMAND TO SAY SOMETHING I THE CONSOLE, AND CALL UPDATE ON ALL LISTENERS TO THIS SHARABLE
        }
        **/
    }
    
    public final void submit(JSONObject data) {
    	
    	// process the data in the JSON object
    	
    }

    public final String create(JSONObject json) {
    
    	System.out.println("[ServerSideAPI.create] "+json);
    	
    	try {
    		
    		String uuid = json.getString("uuid");
    		JSONObject model = json.getJSONObject("model");
    	
    		// here we need to forward the model to the Persistor
    		Persistor persistor = new Persistor(new MongoDBConnection());
    		return persistor.save(model);
    		
    	} catch(Exception e) {
    		e.printStackTrace();
    		return e.getMessage();
    	}
    }
    
    public final JSONObject create(String schemaRef, String jsonArgs) {
    	/***
        try {
            Schema schema = (Schema) resolveToDatum(schemaRef);
            if (schema == null) {
                return null;
            }
            JSONObject json = resolveToJSONObject(jsonArgs);

            //Note:  somehow we need to get the Person from the client
            Person author = Person.getAdmin();

            Instance out = Instance.create(author, schema, json.toString());
            return out.toJSON();
        } catch (ClassCastException e) {
            Logger.log(Logger.Level.WARN, "not a Schema: " + schemaRef, e);
        }
        ***/
        return new JSONObject();
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
    	/**
        try {
            Sharable sharable = (Sharable) resolveToDatum(sharableRef);

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
            Function assistant = (Function) resolveToDatum(asstRef);
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
            Persistor.get().persistDatum(mind);
        } catch (Exception e) {
            Logger.log(Logger.Level.WARN, "", e);
            e.printStackTrace();
            //REVERT THE MIND
            //SEND CLIENT MESSAGE TELLING THAT IT FAILED
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
            View view = (View) resolveToDatum(viewId);
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
            Persistor.get().persistDatum(mind);
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
            Trail trail = (Trail) resolveToDatum(trailRef);
            Person student = this.getPerson();

            Proctor.get().initiateTrail(student, trail, parentDoo);

        } catch (Exception e) {
            Logger.log(Logger.Level.WARN, "", e);
            e.printStackTrace();
        }
        ***/
    }

    public final void gradeQuiz(String answerData) {
    }

    public final void learn(String nativeCmd, String jsCmd) {
        Interpreter.get().learnNative(nativeCmd, jsCmd);
    }

    public final void testout() {
    	/**
        List<Sharable> listy = new ArrayList<Sharable>();
        listy.add(Person.getAdmin());
        listy.add(Person.getSchema());
        try {
            JSONArray data = new JSONArray();
            for (Sharable item : listy) {
                data.put(item.toJSON());
            }

            JSONObject msg = new JSONObject();
            msg.put("collect", data);
            Communicator.get().sendClientMessage(socket_id,
                    SendChannels.commandList,
                    msg.toString());
        } catch (Exception ex) {
            Logger.log(Logger.Level.WARN,
                    "Error sending Sharables to the client");
            ex.printStackTrace();
        }
        **/
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
    public final void callback(String dooID) {
    	/**
        System.out.println("calling back!");
        //Grab the Doo and perhaps perform any commands embedded in there

        //Terminate it from the Hopper
        Hopper.get().terminate(dooID);
        **/
    }

    /***                                     ***\
     ********  Client Command Factories  ********
    \*                                        */
    /**
     * Instruct the client to receive a List of Sharables and log
     * those objects into the clientside collector
     * 
     * @param sharables
     * @return
     * @throws Exception 
     */
    public static JSONObject makeCollect(List<Sharable> sharables) throws Exception {
        JSONObject out = new JSONObject();
        JSONArray data = new JSONArray();
        for (Sharable sharable : sharables) {
            data.put(sharable.toJSON());
        }
        out.put("data", data);
        out.put("command", "collect");

        return out;
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
    public static JSONObject makeShowWidget(Widget widget, String position) throws Exception {
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

    //JCA:  I DON'T KNOW ABOUT THIS ONE...THAT REQUIRES SOME THINKING AS TO WHETHER THIS IS AN OPTION, AND IF SO, HOW IT OCCURS
    //THERE MINIMALLY SHOULD BE SOME DOOS INVOLVED WITH A LOGGING PROCESS
    public void log(String logLevel, String msg) {
        Logger.log(Enum.valueOf(Logger.Level.class, logLevel),
                "@@@@@@@@@@ " + msg);
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
    private Datum resolveToDatum(Object datumRef) {
        Datum out = null;
        Logger.log(Logger.Level.INFO, "resolveToDatum has " + datumRef);

        //If it's a Datum object
        try {
            Datum aschmema = (Datum) datumRef;
            String uuid = aschmema.getId();

            //Fetch the real version of the Datum
            out = Collector.get().getDatum(uuid);
            if (out != null) {
                return out;
            }
        } catch (Exception e) {
            Logger.log(Logger.Level.INFO, "It wasn't a Datum object");
        }

        //If it's a String
        try {
            String str = (String) datumRef;

            //See if it's a uuid String
            out = Collector.get().getDatum(str);
            if (out != null) {
                return out;
            }

            //See if it's a json string representation of a Datum
            JSONTokener tokener = new JSONTokener(str);
            JSONObject jsonobj = new JSONObject(tokener);

            String uuid = jsonobj.getString("uuid");
            out = Collector.get().getDatum(uuid);
            if (out != null) {
                return out;
            }

            //See if it's the proper value of a NameField of some Datum
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
    public static final JSONObject resolveToJSONObject(Object jsonArgs) {
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
            Datum d = Collector.get().getDatum(uuid);
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
            Sharable sharable = (Sharable) Collector.get().getDatum(uuid);
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

    
    /***                       ***\
     ********  VARIABLES  ********
    \*                         */
    /*
    JCA:  I SORTOVE DON'T LIKE THIS CALLBACK BUSINESS ANYMORE.  THERE SHOULD BE A DOO IN A DOO
    HOPPER BEHAVING AS THE CALLBACK, NOT A GENERIC CALLBACK.  THIS DOESN'T ALLOW A
    LOGGING OF THE CHAIN OF EVENTS THAT LED TO THE CALLBACK.
    
    ALSO, EVERYTHING GOING ON WITH EXECUTION IS GOING TO BE SINGLE THREADED WITH THE NON-BLOCKING LOOP
    BUSINESS HANDLING ONLY INCOMING MESSAGES.  SO, THIS IS GETTING EXECUTED IN A FULLY SYNCHRONIZED
    AND UNINTERUPTED MANNER THROUGH TO COMPLETION WHATEVER IT IS.  SO, IT'S NOT GOING TO BE AN
    ASYNCHRONOUS PATTERN
    
    WELL, OK, I DON'T KNOW.  MAYBE YOU DO WANT EXECUTION TO BE INTERNALLY ASYNCHRONOUS TOO, OR PERHAPS
    YOU WANT ASSISTANTS TO BE MULTITHREADED SO THAT YOU CAN ISSUE A SET OF COMMANDS THAT EACH INVOLVE
    SAY WAITING FOR A PUBMED RESPONSE.  THAT WOULD BE INTERNALLY ASYNCHRONOUS.  I'M NOT SURE, SO I'M
    GOING TO LEAVE THIS FOR SOMEONE ELSE TO RESOLVE.  YOU DON'T WANT THE ENTIRE CLOTHO HELD UP BECAUSE
    IT'S WAITING FOR AN ATTEMPT TO CONNECT TO ANOTHER SERVER TO TIME-OUT, WHICH YOU'D GET IF IT WAS
    ENTIRELY SINGLE-THREADED.  IT PROBABLY DOES ALL HAVE TO BE ASYNCHRONOUS, BUT IT SHOULDN'T BE VANILLA
    CALLBACKS, IT SHOULD BE DOOS.  THE ABILITY TO TRACK AND NOT LOSE DO'S ULTIMATELY WILL BE VERY IMPORTANT
    FOR ENABLING LOGGING, MAINTAINING STATISTICS ABOUT THINGS, TRACKING HISTORY OF EVENTS FOR SECURITY AND
    LEGAL REASONS.  DOOS ARE FOR ALL SORTS OF STUFF, SO THEY NEED TO BE PASSED AROUND CONSTANTLY.*/
    private Callback callback = new Callback() {

        @Override
        public void onSuccess(JSONObject outputData) {
            Logger.log(Logger.Level.INFO,
                    "All done! outputData = " + String.valueOf(outputData));
        }

        @Override
        public void onFailure(Throwable e) {
            Logger.log(Logger.Level.WARN, "All done! but failed", e);
        }
    };

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
    
    //private transient Person person;
    private static final AutoComplete completer = new AutoComplete();
}

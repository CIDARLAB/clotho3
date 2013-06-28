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

package org.clothocad.core.layers.communication.mind;

import java.io.File;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;

import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;

import org.clothocad.core.aspects.Aspect;
import org.clothocad.core.datums.Doo;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.layers.communication.ServerSideAPI;
import org.clothocad.core.layers.communication.connection.ClientConnection;
import org.clothocad.core.util.FileUtils;
import org.clothocad.model.Person;
import org.json.JSONException;
import org.json.JSONObject;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * A 'Mind' is the thing on the server that connects up a User and a Client.  Minds talk
 * to the Communicator (the singleton, to it's many instances).
 *
 * A User might use multiple clients and do different things on different machines
 *
 * A Client might be used by a different humans, and one human might be using
 * several Persons
 *
 * So, a Mind is where things get placed as being specific to one Person and one Client
 *
 * All communication between a client and a user goes through a Mind (it will
 * contextualize the session using its namespace, so this uuid is what goes into queries as a hint).
 * 
 * It also maintains the physical configuration of the client and persists that.
 *
 * @author John Christopher Anderson
 */

public final class Mind 
		extends ObjBase
		implements Aspect {
	
    static final Logger logger = LoggerFactory.getLogger(Mind.class);
    
    public static Mind create(ClientConnection connection) {
        Mind mind = new Mind();
        mind.connection = connection;        
        mind.ssAPI = new ServerSideAPI();
        mind.ssAPI.setMind(mind);
        return mind;
    }
    /**
     * This is not a serverside API method
     * @param client
     * @param person
     * @return 
     */
    public static Mind create(Person person) {
        /* TODO check if this Mind already exists */
        //Create the new Mind object
        Mind out = new Mind();
        out.personId = person.getId();
        //out.save();
        return out;
    }

    //THIS SHOULD BE CHANGED TO PRIVATE, BUT TESTERS NEED IT FOR NOW
    public Mind() {
        config = new ClientConfiguration();
        engine = getEngine();
    }

    /**
     * Serverside API relayed method to 'clear the mind'
     */
    public final void clear() {
        engine = null;
        getEngine();
    }

    /* Similar to runCommand but does not "defuzzify".
     * Used for the "serverEval" channel. */
    public synchronized boolean eval(String socket_id, String cmd) {
        try {
            getEngine().eval(cmd);
            return true;
        } catch (ScriptException e) {
            return false;
        }
    }

    /**
     * Try to run the command as if it were a proper script.  If fails,
     * try it again after injecting some objects that might not be in 
     * the engine yet.  If all that fails, return, then abort and disambiguate.
     * @param socket_id
     * @param cmd
     * @return 
     */
    /* TODO: race condition. ScriptEngine execution needs to be serialized. */
    public synchronized boolean runCommand(String cmd) {
        try {
            getEngine().eval(cmd);
        } catch (ScriptException e) {
            return false;
        }
        
        return true;
    }

    /* client has new tab--update config */
    public void linkPage(String socket_id,
                         String ephemeral_link_page_id,
                         PageMode page_mode) {
        config.addPage(socket_id, page_mode);
        //XXX: this.save();
        if (doos.containsKey(ephemeral_link_page_id)) {
            AddPageDoo add_page_doo =
                (AddPageDoo) doos.get(ephemeral_link_page_id);
            if (add_page_doo.getPageMode().compareTo(page_mode) == 0) {
                add_page_doo.terminate();
                /* TODO: add_page_doo is done
                 * call some sort of callback for the caller of addPage
                 */
            } else {
                logger.error(
                           "linkNewPage got page mode mismatch: {}, {}",
                           page_mode.name(),
                           add_page_doo.getPageMode().name());
            }
        }
    }
    
    /* Instruct client to close a tab */
    public void removePage(String socket_id) {
         /* TODO: DO WHATEVER IT TAKES TO programmatically remove THE TAB */
    }

    /* client lost a page--update config */
    public void unlinkPage(String socket_id) {
        config.removePage(socket_id);
        //XXX: this.save();
    }
    
    public void setVisible(Page page, boolean new_visibility) {
        /* TODO: DO WHATEVER IT TAKES TO programmatically show or hide the tab */
        
        page.toggleVisible(new_visibility ^ page.toggleVisible(false));
        //XXX: this.save();
    }

    public Iterable<Map<String, String>> getPageSummary() {
        List<Map<String, String>> out = new ArrayList<Map<String, String>>();
        for (String socket_id : config.getSocketIDs()) {
            Map<String, String> out_item = new HashMap<String, String>();
            out_item.put("socket_id", socket_id);
            out_item.put("name", "(No Name)");
            out_item.put("mode",
                         config.getPage(socket_id).getPageMode().name());
            out.add(out_item);
        }
        return out;
    }    

    public void removeWidget(Page page, Widget widget) {
        /* TODO: remove the widget from the client's page */
        page.removeWidget(widget.getId());
        //XXX: this.save();
    }

    public void unlinkWidget(String socket_id, String widget_id) {
        Page target_page = config.getPage(socket_id);
        if (target_page == null) {
            logger.error("can't get page with socket_id: " + socket_id);
            return;
        }
        target_page.removeWidget(widget_id);
        //XXX: this.save();
    }

    public void moveWidget(Widget widget,
                           Page page,
                           Double posx,
                           Double posy,
                           Double width,
                           Double height) {
        //DO WHATEVER IT TAKES TO programmatically move the widget inside the client
        
        widget.setPositionAbsolute(posx, posy);
        widget.setDimensions(width, Widget.SizeType.RELATIVE, height, Widget.SizeType.RELATIVE);
        
        //XXX: this.save();
    }

    /* TODO: this shouldn't be exposed */
    ClientConfiguration getClientConfiguration() {
        return config;
    }

    /* TODO: this shouldn't be exposed */
    ScriptEngine getMyEngine() {
        return engine;
    }

    public String getPersonId() {
        return personId;
    }
    
    /**
     * THIS ALL FEELS A LITTLE FISHY--IT PROBABLY NEEDS TO SAVE THIS DATA CHANGE
     * NEED TO FIGURE OUT WHY THERE IS A NEED TO SET THIS IN NameSpace
     * @param newid 
     */
    public void setPersonId(String newid) {
        this.personId = newid;
    }
    
    /**
     * Return the scripting engine
     * (or a new instance if it doesn't already exist)
     *
     * Injects a new ServerScriptAPI with `socket_id`
     * `socket_id` may be null if this request is not bound to a page.
     */
    public ScriptEngine getEngine() {
        if (engine == null) {
            engine = new ScriptEngineManager().getEngineByName("JavaScript");
            ServerSideAPI api = getAPI();
            engine.put("clothoJava", api);
            
            try {
                engine.eval(initializationScript);
            } catch (ScriptException ex) {
                System.out.println("Error running initialization Script!");
                ex.printStackTrace();
            }
        }
        return engine;
    }

    public ServerSideAPI getAPI() {
        if (ssAPI == null) {
            ssAPI = new ServerSideAPI(this);
        }
        return ssAPI;
    }
    
    public Person getPerson() {
    	return this.person;
    }
    
    public void setClientConnection(ClientConnection conn) {
        connection = conn;
    }
    
    public List<JSONObject> getLastCommands() {
        return this.lastCommands;
    }
    
    public void addLastCommand(String str) {
        try {
            JSONObject obj = new JSONObject();
            obj.put("text", str);
            obj.put("type", "phrase");
            lastCommands.add(obj);
        } catch (JSONException ex) {
            ex.printStackTrace();
        }
    }
    
    public void addLastCommand(JSONObject json) {
        lastCommands.add(json);
    }
    
    private Person person;
    private String personId;
    private transient List<JSONObject> lastCommands = new ArrayList<JSONObject>();
    private transient List<String> lastSharables = new ArrayList<String>();
    private transient ScriptEngine engine;
    private transient ServerSideAPI ssAPI;
    private transient ClientConnection connection;
    
    private Map<Date, String> commandHistory;
    
    private static final int MAX_SIZE_LAST_SHARABLES = 50;

    /* How the client is currently set up */
    private ClientConfiguration config;

    /* Whether it should persist the visual configuration */
    private boolean useInitialConfiguration = false;

    /* Maps "Doo ID" to Doo object */
    
    //JCA:  THIS IS PROBLEMATIC, IT IS HOLDING STATE FOR A DOO, NOT LOOSE-COUPLED
    //NOT SURE WHY THIS WOULD BE HERE INSTEAD OF THE HOPPER, I THINK THIS IS PROBABLY WRONG
    private Map<String, Doo> doos = new HashMap<String, Doo>();

    private String id;

public void SUPERILLEGAL_SETUUID(String string) {
    id = string;
}

    public ClientConfiguration getConfig() {
        return config;
    }

    public ClientConnection getClientConnection() {
       return this.connection;
   }

    public String pullUUIDFromNamespace(String str) {
        System.out.println("JCA:  Need to implement mapping names/namespace tokens to UUIDs of Sharables");
        return null;
    }
    
    private static  String initializationScript;
    static {
        System.out.println("Someone:  Probably should pull this from elsewhere");
        initializationScript = FileUtils.readFile("js_engine_initiation.js");
    }

    public List<String> getRecentSharables() {
        return this.lastSharables;
    }

    public void addLastSharable(String lastid) {
        if(this.lastSharables.contains(lastid)) {
            lastSharables.remove(lastid);
        }
        lastSharables.add(lastid);
        if(lastSharables.size() > MAX_SIZE_LAST_SHARABLES) {
            lastSharables.remove(MAX_SIZE_LAST_SHARABLES);
        }
    }

}

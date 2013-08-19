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

import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.script.ScriptException;
import lombok.Getter;
import lombok.Setter;

import org.clothocad.core.aspects.Aspect;
import org.clothocad.core.datums.Doo;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.layers.communication.Channel;
import org.clothocad.core.layers.communication.Message;
import org.clothocad.core.layers.communication.ScriptAPI;
import org.clothocad.core.layers.communication.connection.ClientConnection;
import org.clothocad.core.layers.execution.MetaEngine;
import org.clothocad.core.util.FileUtils;
import org.clothocad.model.Person;
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
        out.personId = person.getUUID().toString();
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
    public synchronized Object eval(String cmd, ScriptAPI api) throws ScriptException {
            return getEngine().eval(cmd, Language.JAVASCRIPT, api);    
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
    public synchronized Object runCommand(String cmd, ScriptAPI api) throws ScriptException {
        return eval(cmd, api);
//                runCommandWithNamespace(cmd);
    }

    /* client has new tab--update config */
    public void linkPage(String socket_id,
                         String ephemeral_link_page_id,
                         PageMode page_mode) {
        throw new UnsupportedOperationException();
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
        List<Map<String, String>> out = new ArrayList<>();
        for (String socket_id : config.getSocketIDs()) {
            Map<String, String> out_item = new HashMap<>();
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
     */
    public MetaEngine getEngine() {
        if (engine == null) {
            engine = new MetaEngine();
        }
        return engine;
    }
    
    public Person getPerson() {
    	return this.person;
    }
    
    public List<Message> getLastCommands() {
        return this.lastCommands;
    }

    public void addLastCommand(Message command){
        lastCommands.add(command);
    }
    
    public void addLastCommand(Channel channel, Object data){
        lastCommands.add(new Message(channel, data));
    }
    
    private Person person;
    private String personId;
    private transient List<Message> lastCommands = new ArrayList<>();
    private transient MetaEngine engine;
    private transient List<String> lastSharables = new ArrayList<>();
    @Getter @Setter
    private transient ClientConnection connection;
    
    private Map<Date, String> commandHistory;
    
    private static final int MAX_SIZE_LAST_SHARABLES = 50;

    /* How the client is currently set up */
    @Getter
    private ClientConfiguration config;

    /* Whether it should persist the visual configuration */
    private boolean useInitialConfiguration = false;

    /* Maps "Doo ID" to Doo object */
    
    //JCA:  THIS IS PROBLEMATIC, IT IS HOLDING STATE FOR A DOO, NOT LOOSE-COUPLED
    //NOT SURE WHY THIS WOULD BE HERE INSTEAD OF THE HOPPER, I THINK THIS IS PROBABLY WRONG
    private Map<String, Doo> doos = new HashMap<>();

    private String id;

public void SUPERILLEGAL_SETUUID(String string) {
    id = string;
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
    
    public Object evalFunction(String code, String name, List args, ScriptAPI api) throws ScriptException {
        return getEngine().invoke(code, name, args);
    }

}

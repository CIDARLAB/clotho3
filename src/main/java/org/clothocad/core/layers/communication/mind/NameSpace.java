package org.clothocad.core.layers.communication.mind;

import java.util.Map;
import java.util.HashMap;
import java.util.TreeMap;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/* TODO: Move this functionality (and the related methods in Mind)
 *       to Interpreter.
 *
 * Maintains the list of "var"-token style commands issued into the
 * script engine over your entire history.  It has a fixed size of
 * 10000 items, but it will drop old tokens over time to let in new ones.
 */
class NameSpace {
    static final Logger logger = LoggerFactory.getLogger(NameSpace.class);
    
    public void put(String token, String cmd) {
        //If the token bank hasn't seen this token
        if(!tokenBank.containsKey(token)) {
            //If it hasn't reached its maximum size
            if(tokenBank.size()<10000) {
                tokenBank.put(token, cmd);
                int lasttime = 0;
                if(!timeline.isEmpty()) {
                    lasttime = timeline.lastKey();
                }
                
                //Renumber indexes if they get too long
                if(lasttime>99000) {
                    renumber();
                    put(token, cmd);
                    return;
                }
                
                timeline.put(lasttime + 1, token);
            } 
            
            //If it has reached its maximum size, remove oldest and try again
            else {
                int firsttime = timeline.firstKey();
                String oldtoken = timeline.get(firsttime);
                tokenBank.remove(oldtoken);
                timeline.remove(firsttime);
                put(token, cmd);
            }
        } 
        
        //If it has seen the token, update that
        else {
            tokenBank.put(token, cmd);
        }
    }
    
    public boolean containsKey(String word) {
        return tokenBank.containsKey(word);
    }

    public String get(String word) {
        return tokenBank.get(word);
    }

    public void learn(String token, String cmd, Mind mind) {
        Mind newMind = new Mind();  //Won't be persisted
        newMind.setPersonId(mind.getPersonId());
        
        //Issue the command in a clean engine; if it fails, don't try to learn
        try {
            /* Not sure if getting a new engine will be different from the
             * already instantiate engine */
//            newMind.getMyEngine().eval(cmd);  //this method is problematic as it is now taking in arguments for show and trying to show them
            put(token, cmd);
            //XXX: mind.save();
        } catch (Exception e) {
            logger.warn( "won't learn command: {}");
            return;
        }
    }

    private void renumber() {
        TreeMap<Integer, String>  out = new TreeMap<Integer, String>();
        int count = 0;
        for(int key : timeline.navigableKeySet()) {
            out.put(count, timeline.get(key));
            count++;
        }
        timeline = out;
    }

    Map<String, String> tokenBank = new HashMap<String, String>();
    TreeMap<Integer, String> timeline = new TreeMap<Integer, String>();
}

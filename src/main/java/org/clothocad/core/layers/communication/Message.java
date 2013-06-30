/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.layers.communication;

import java.util.List;
import java.util.Map;

/**
 *
 * @author spaige
 */
public class Message {
    
    public Message(Channel channel, Object data){
        this(channel, data, null);
    }
    
    public Message(Channel channel, Object data, String id){
        this.channel = channel;
        this.data = data;
        this.requestId = id;
    }
    
    public final Channel channel;
    public final Object data;
    public final String requestId;

    
    //assumes String describes JSONObject with keys "channel" and "data"
    public Message(String messageString) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    
    //XXX: not idempotent if nested lists
    public Message unwrapData(){
        if (data instanceof List && ((List) data).size() == 1)
            return new Message(channel, ((List) data).get(0), requestId);
        return this;
    }
}

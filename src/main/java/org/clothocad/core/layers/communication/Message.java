/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.layers.communication;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import java.util.List;
import java.util.Map;
import org.clothocad.core.util.JSON;

/**
 *
 * @author spaige
 */
public class Message {
    public Message(Channel channel, Object data){
        this(channel, data, null);
    }
    
    @JsonCreator
    public Message(@JsonProperty("channel") Channel channel, 
    @JsonProperty("data") Object data, @JsonProperty("requestId") String requestId){
        this.channel = channel;
        this.data = data;
        this.requestId = requestId;
    }
    
    public final Channel channel;
    public final Object data;
    public final String requestId;

    
    /*//assumes String describes JSONObject with keys "channel","data", and "requestId"
    public Message(String messageString) {
        Map<String,Object> map = JSON.mappify(messageString);
        channel = Channel.valueOf(map.get("channel").toString());
        data = map.get("data");
        requestId = map.get("requestId").toString();
    }*/
    
    
    //XXX: not idempotent if nested lists
    public Message unwrapData(){
        if (data instanceof List && ((List) data).size() == 1)
            return new Message(channel, ((List) data).get(0), requestId);
        return this;
    }
}

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.layers.communication;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import java.util.Map;
import lombok.NoArgsConstructor;

/**
 *
 * @author spaige
 */
//TODO: figure out how to make immutable again w/ Morphia (subclass objectfactory)
@NoArgsConstructor
public class Message {
    public Message(Channel channel, Object data){
        this(channel, data, null);
    }
    

    public Message(@JsonProperty("channel") Channel channel, 
    @JsonProperty("data") Object data, @JsonProperty("requestId") String requestId){
        this.channel = channel;
        this.data = data;
        this.requestId = requestId;
    }

    @JsonCreator
    public Message(Map<String,Object> props){
        this.channel = Channel.valueOf(props.get("channel").toString());
        this.data = props.get("data");
        this.options = (Map) props.get("options");
        this.requestId = props.containsKey("requestId")? props.get("requestId").toString() : null;
    }
    
    public Channel channel;
    public Object data;
    public Map<String,String> options;
    public String requestId;

    
    /*//assumes String describes JSONObject with keys "channel","data", and "requestId"
    public Message(String messageString) {
        Map<String,Object> map = JSON.mappify(messageString);
        channel = Channel.valueOf(map.get("channel").toString());
        data = map.get("data");
        requestId = map.get("requestId").toString();
    }*/   
}

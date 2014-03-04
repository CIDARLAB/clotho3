/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.communication;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.databind.JsonMappingException;
import java.io.IOException;
import java.io.StringWriter;
import java.util.Map;
import lombok.extern.slf4j.Slf4j;
import org.clothocad.core.util.JSON;

/**
 *
 * @author spaige
 */
@Slf4j
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
}

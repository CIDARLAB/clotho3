package org.clothocad.core.communication;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.databind.JsonMappingException;
import java.io.IOException;
import java.io.StringWriter;
import org.clothocad.core.util.JSON;

/**
 * @author spaige
 */
//TODO: figure out how to make immutable again w/ Morphia (subclass objectfactory)
public class Message {
    private final Channel channel;
    private Object data;
    private final String requestId;

    @JsonProperty("channel")
    public Channel
    getChannel() {
        return channel;
    }

    @JsonProperty("data")
    public Object
    getData() {
        return data;
    }

    @JsonProperty("requestId")
    public String
    getRequestId() {
        return requestId;
    }

    public Message(@JsonProperty("channel") Channel channel, 
                   @JsonProperty("data") Object data,
                   @JsonProperty("requestId") String requestId) {
        this.channel = channel;
        this.data = data;
        this.requestId = requestId;
    }

    public String serialize() {
        try {
            StringWriter writer = new StringWriter();
            JSON.mapper.writeValue(writer, this);
            return writer.toString();
        } catch(JsonMappingException ex) {
            //we can safely assume that data is the unserializable component
            data = data.toString();
            return JSON.serialize(this);
        } catch (IOException ex) {
            ex.printStackTrace();
            return ex.getMessage();
        }
    }
}

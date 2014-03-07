package org.clothocad.core.communication;

import com.fasterxml.jackson.annotation.JsonProperty;

/**
 * @author spaige
 */
//TODO: figure out how to make immutable again w/ Morphia (subclass objectfactory)
public class Message {
    private final Channel channel;
    private final Object data;
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
}

package org.clothocad.core.communication;

import com.fasterxml.jackson.annotation.JsonProperty;
import java.util.Map;
import org.clothocad.core.util.MapView;

/**
 * @author spaige
 */
public class Message {
    private final Channel channel;
    private final Object data;
    private final String requestId;
    private final Map<String, String> options;

    @JsonProperty("channel")
    public Channel getChannel() { return channel; }

    @JsonProperty("data")
    public Object getData() { return data; }

    @JsonProperty("requestId")
    public String getRequestId() { return requestId; }

    @JsonProperty("options")
    public Map<String, String> getOptions() { return options; }

    public Message(@JsonProperty("channel") Channel channel,
                   @JsonProperty("data") Object data,
                   @JsonProperty("requestId") String requestId,
                   @JsonProperty("options") Map<String, String> options) {
        this.channel = channel;
        this.data = data;
        this.requestId = requestId;
        this.options = MapView.wrap(options);
    }
}

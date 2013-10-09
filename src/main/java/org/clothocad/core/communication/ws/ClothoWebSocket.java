package org.clothocad.core.communication.ws;

import com.fasterxml.jackson.core.JsonParseException;
import com.fasterxml.jackson.databind.JsonMappingException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.Callable;
import lombok.extern.slf4j.Slf4j;
import org.apache.shiro.SecurityUtils;
import org.apache.shiro.subject.Subject;
import org.clothocad.core.communication.Channel;
import org.clothocad.core.communication.Message;
import org.clothocad.core.communication.Router;
import org.clothocad.core.communication.ClientConnection;
import org.clothocad.core.util.JSON;
import org.eclipse.jetty.websocket.WebSocket;

@Slf4j
public class ClothoWebSocket 
	extends ClientConnection 
	implements WebSocket.OnTextMessage {

	private WebSocket.Connection connection;
        
            private Subject subject;
    private final Router router;

    private class CallRouter implements Callable {
        private final Message message;
        private final ClientConnection connection;
        CallRouter(ClientConnection connection, Message message){
            this.message = message;
            this.connection = connection;
        }

        @Override
        public Object call() throws Exception {
            router.receiveMessage(connection, message);
            return null;
        }
        
    }
	
	public ClothoWebSocket(String id, Router router) {
		super(id);
                this.router = router;
	}
	
	@Override
	public void onClose(int closeCode, String message) {
	}

    @Override
    public void send(Message msg) {
        try {
            connection.sendMessage(msg.serialize());
            log.trace("sent: {}", msg.serialize());
        } catch (IOException ex) {
            log.error("Cannot send message", ex);
        }
    }

    @Override
    public void onMessage(String messageString) {
        log.trace("Websocket #{} recieved message {}", this.getId(), messageString);
        try {
            Message message = JSON.mapper.readValue(messageString, Message.class);
            subject.execute(new CallRouter(this, message));
        } catch (JsonParseException ex) {
            log.error("Websocket #{} recived malformed message: {}", this.getId(), messageString);
        } catch (JsonMappingException ex) {
            throw new RuntimeException(ex);
        } catch (IOException ex) {
        }   catch (Exception ex) {
                throw new RuntimeException(ex);
        }
    }

	public boolean isOpen() {
		return connection.isOpen();
	}

	@Override
	public void onOpen(Connection connection) {
		this.connection = connection;
                //Close out after 1 hour idle time
                connection.setMaxIdleTime(3600000); 

		subject = SecurityUtils.getSubject();
		
	}
        
    
    @Override
    public void deregister(Channel channel, String requestId) {
        Map<String, String> message = new HashMap<>();
        message.put("channel", channel.toString());
        message.put("requestId", requestId);
        try {
            connection.sendMessage(JSON.serialize(message));
            log.trace("sent deregister: {}", JSON.serialize(message));
        } catch (IOException ex) {
            log.error("cannot send deregister message:{}", ex.getMessage());
        }
    }
}

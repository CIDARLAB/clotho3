package org.clothocad.core.layers.communication.connection.ws;

import com.fasterxml.jackson.core.JsonParseException;
import com.fasterxml.jackson.databind.JsonMappingException;
import java.io.IOException;
import lombok.extern.slf4j.Slf4j;

import org.clothocad.core.layers.communication.Message;
import org.clothocad.core.layers.communication.Router;
import org.clothocad.core.layers.communication.connection.ClientConnection;
import org.clothocad.core.util.JSON;
import org.eclipse.jetty.websocket.WebSocket;

@Slf4j
public class ClothoWebSocket 
	extends ClientConnection 
	implements WebSocket.OnTextMessage {

	private Connection connection;
	
	public ClothoWebSocket(String id) {
		super(id);
	}
	
	@Override
	public void onClose(int closeCode, String message) {
            
	}

    @Override
    public void send(Message msg) {
            try {
                connection.sendMessage(JSON.serialize(msg));
            } catch (IOException ex) {
                log.error("Cannot send message", ex);
            }
    }
        
        

	@Override
	public void onMessage(String messageString) {
            log.trace("Websocket #{} recieved message {}", this.getId(), messageString);
                try {
                    Message message = JSON.mapper.readValue(messageString, Message.class);
                    Router.get().receiveMessage(this, message);
                }  catch (JsonParseException ex) {
                    log.error("Websocket #{} recived malformed message: {}", this.getId(), messageString);
                } catch (JsonMappingException ex) {
                    throw new RuntimeException(ex);
                } catch (IOException ex) {
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

		//WebSocketTable.put(this.getId(), this);
		
		// TODO: store the connection information into the Mind
		
	}
}

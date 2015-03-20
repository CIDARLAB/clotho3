package org.clothocad.core.communication.ws;

import com.fasterxml.jackson.core.JsonParseException;
import com.fasterxml.jackson.databind.JsonMappingException;
import java.io.IOException;
import java.net.URI;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import lombok.extern.slf4j.Slf4j;
import org.apache.shiro.SecurityUtils;
import org.apache.shiro.subject.Subject;
import org.clothocad.core.communication.Channel;
import org.clothocad.core.communication.Message;
import org.clothocad.core.communication.Router;
import org.clothocad.core.communication.ClientConnection;
import org.clothocad.core.util.JSON;
import org.eclipse.jetty.websocket.WebSocket;
import org.eclipse.jetty.websocket.WebSocketClient;
import org.eclipse.jetty.websocket.WebSocketClientFactory;

@Slf4j
public class ClothoWebSocketServer
        extends ClientConnection
        implements WebSocket.OnTextMessage {

    private WebSocket.Connection connection;
    private Subject subject;
    private final Router router;
    
    
    private HashMap<String,WebSocket.Connection > serverConnections;    //uri to connection
    private HashMap<String, Message> serverResponses;    //uri to message
    private HashMap<String, Boolean> serverResponseReceived;    //uri to if message was received
    private String uri;         //uri to create the connection
    
    private class CallRouter implements Callable {

        private final Message message;
        private final ClientConnection connection;

        CallRouter(ClientConnection connection, Message message) {
            this.message = message;
            this.connection = connection;
        }

        @Override
        public Object call() throws Exception {
            router.receiveMessage(connection, message);
            return null;
        }
    }
    private static ClothoWebSocketServer clothoWebSocket; //by Jin
    
    private ClothoWebSocketServer(String id, Router router) { //changed to private
        super(id);
        this.router = router;
        
        serverConnections = new HashMap<String,WebSocket.Connection >();    //uri to connection
        serverResponses = new HashMap<String, Message>();    //uri to message
        serverResponseReceived = new HashMap<String, Boolean>();    //uri to if message was received
    }
    public static ClothoWebSocketServer getInstance(){
        return clothoWebSocket;
    }
    //method to get the instance of clothowebsocket
    public static ClothoWebSocketServer getInstance(String id, Router router){
        if(clothoWebSocket == null){
            clothoWebSocket = new ClothoWebSocketServer(id, router);
        }
        return clothoWebSocket;
    }
    @Override
    public void onClose(int closeCode, String message) {
    }

    @Override
    public void send(Message msg) {
        try {
            String messageString = JSON.serializeForExternal(msg);
            connection.sendMessage(messageString);
            log.trace("sent: {}", messageString);
        } catch (IOException ex) {
            log.error("Cannot send message", ex);
        }
    }
    public WebSocket.Connection getConnection(String uri){
        this.uri = uri;
        if (serverConnections.containsKey(uri)){
            return serverConnections.get(uri);
        }else{
            Connection newCon = createConnection(uri);
            serverConnections.put(uri, newCon);
            serverResponseReceived.put(uri, Boolean.FALSE);
            return newCon;
        }
    }
    
    @Override
    public void onMessage(String messageString) {
        log.trace("Websocket #{} recieved message {}", this.getId(), messageString);
        try {
            Message message = JSON.mapper.readValue(messageString, Message.class);

            if(message.getOptions() == null && messageString.contains("options")){  //terminating condition
                System.out.println("Server received a response, not a message from the client");
                serverResponses.put(uri, message);
                serverResponseReceived.put(uri, Boolean.TRUE);
            }else{  
                System.out.println("Server received a query");
                subject.execute(new CallRouter(this, message));
            }
        } catch (JsonParseException ex) {
            log.error("Websocket #{} recived malformed message: {}", this.getId(), messageString);
        } catch (JsonMappingException ex) {
            throw new RuntimeException(ex);
        } catch (IOException ex) {
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
    }

    public boolean isOpen() {
        return connection.isOpen();
    }

    @Override
    public void onOpen(Connection connection) {
        log.debug("New connection opened. Connection id is {}.", this.getId());
        // before storing, send a message back to see what kind of connection it is
        
        this.connection = connection;   //call router automatically makes the connection

        connection.setMaxIdleTime(3600000);

        subject = SecurityUtils.getSubject();
    }


    @Override
    public void deregister(Channel channel, String requestId) {
        String messageString = String.format("{\"channel\": \"%s\", \"requestId\": \"%s\"}", channel, requestId);
        try {
            connection.sendMessage(messageString);
            log.trace("sent deregister: {}", messageString);
        } catch (IOException ex) {
            log.error("cannot send deregister message:{}", ex.getMessage());
        }
    }
   
    private Connection createConnection(String destUri){
        try{
            WebSocketClientFactory factory = new WebSocketClientFactory();
            factory.start();
            WebSocketClient wsClient = factory.newWebSocketClient();
            URI uri = new URI(destUri);
            Future fut = wsClient.open(uri, clothoWebSocket); 
            return (Connection) fut.get(10, TimeUnit.SECONDS);
        }catch(Throwable t){
            t.printStackTrace();
        }
        return null;
    }
    public boolean gotMessage(String uri){
        return serverResponseReceived.get(uri);
    }
    public Message getMessage(String uri){
        return serverResponses.get(uri);
    }
}

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.webserver.jetty;

import java.io.IOException;
import org.clothocad.core.communication.Router;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Component;
import org.springframework.web.socket.CloseStatus;
import org.springframework.web.socket.WebSocketHandler;
import org.springframework.web.socket.WebSocketMessage;
import org.springframework.web.socket.WebSocketSession;
import org.springframework.web.socket.adapter.jetty.JettyWebSocketSession;

/**
 *
 * @author david
 */
@Component
public class ClothoWSHandler implements WebSocketHandler {

    /*
    
    JettyWebSocketHandlerAdapter exists, and so does JettyWebSocketSession. They are websocket adapters for Spring => Jetty 9
    
    We're currently trying to figure out how we can hand over the Spring Websocket & sessions to Jetty using the adapter library.
    
    
    Very little documentation or samples on this library and implementing it. 
    JettyWebSocketHandler Adapter does not implement WebSocket Handler, so we have to register the Spring websocket.
    We also don't know if WebSocketSession that gets passed to the Spring WebSocketHandler is even a JettyWebSocketSession, but yolo we casted it.
    
    
    Able to handle the message up to the point where the router actually receives the message.
    
     */
    private Router router;
    private ClothoJettyHandler jettyHandler;
    private JettyWebSocketSession jettySession;
    

    private static Logger logger = LoggerFactory.getLogger(ClothoWSHandler.class);

    public void setRouter(Router router) {
        this.router = router;
    }

    @Override
    public void afterConnectionEstablished(WebSocketSession session) {
        try {

            //SpringClothoStarter adapts WebSocket connections to Jetty 9 API when registering handlers with the Jetty handshaker
            //We can cast spring WebSocketSessions to JettyWebSocketSessions, and initialize it with itself.
            jettySession = (JettyWebSocketSession) session;
            jettySession.initializeNativeSession(jettySession.getNativeSession());
            
            //Setup our custom handler and populate the router, which was magically created by Spring using our annotated classes in Router and its dependencies.       
            jettyHandler = new ClothoJettyHandler(this, jettySession);
            jettyHandler.setRouter(router);
            
            //Let the jetty web socket commence!
            jettyHandler.onOpen(jettySession.getNativeSession());
        } catch (Throwable ex) {
            logger.debug(this.getClass().getName() + ": afterConnectionEstablishedException - " + ex.toString() + ": " + ex.getMessage() + ": " + ex.getLocalizedMessage());
        
        }
    }

    @Override
    public void handleMessage(WebSocketSession session, WebSocketMessage<?> message) throws IOException {
        try {
            logger.debug("Handler open: " + jettyHandler.isOpen());
            logger.debug("Session still open: " + jettySession.isOpen());
            logger.debug("Message being sent to router: " + message.getPayload().toString());
            jettyHandler.onMessage(message.getPayload().toString());
        } catch (Throwable ex) {
            logger.debug(this.getClass().getName() + ": handleMessageException - " + ex.getCause().toString()+ ": " + ex.getMessage() + ": " + ex.getLocalizedMessage());
            session.close();
        }
    }

    @Override
    public void handleTransportError(WebSocketSession session, Throwable exception
    ) {

    }

    @Override
    public void afterConnectionClosed(WebSocketSession session, CloseStatus closeStatus) throws IOException {
        logger.debug("Closing the session: " + session.getId());
        jettyHandler.onClose(closeStatus.getCode(), "Closing Connection");
        session.close(closeStatus);
    }

    @Override
    public boolean supportsPartialMessages() {
        return true;
    }

}

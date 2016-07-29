/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.webserver.jetty;

import java.util.Map;
import org.clothocad.core.communication.Channel;
import org.clothocad.core.communication.Message;
import org.clothocad.core.communication.MessageOption;
import org.clothocad.core.communication.Router;
import org.clothocad.core.communication.ServerSideAPI;
import org.clothocad.core.communication.ws.ClothoWebSocket;
import org.clothocad.core.util.JSON;
import org.eclipse.jetty.websocket.api.Session;
import org.json.JSONObject;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;
import org.springframework.web.socket.CloseStatus;
import org.springframework.web.socket.TextMessage;
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


    private static Logger logger = LoggerFactory.getLogger(ClothoWSHandler.class);

    public void setRouter(Router router) {
        this.router = router;
    }

    @Override
    public void afterConnectionEstablished(WebSocketSession session) throws Exception {
        
        jettyHandler = new ClothoJettyHandler(this, (JettyWebSocketSession) session);
        
        jettyHandler.setRouter(router);
    }

    @Override
    public void handleMessage(WebSocketSession session, WebSocketMessage<?> message) throws Exception {
        
    }

    @Override
    public void handleTransportError(WebSocketSession session, Throwable exception) throws Exception {

    }

    @Override
    public void afterConnectionClosed(WebSocketSession session, CloseStatus closeStatus) throws Exception {
        logger.debug("Closing the session: " + session.getId());
        session.close(closeStatus);
    }

    @Override
    public boolean supportsPartialMessages() {
        return true;
    }

}

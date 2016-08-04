/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.webserver.jetty;

import com.fasterxml.jackson.core.JsonParseException;
import com.fasterxml.jackson.databind.JsonMappingException;
import java.io.IOException;
import java.util.concurrent.Callable;
import org.apache.shiro.SecurityUtils;
import org.apache.shiro.subject.Subject;
import org.clothocad.core.communication.ClientConnection;
import org.clothocad.core.communication.Message;
import org.clothocad.core.communication.Router;
import org.clothocad.core.communication.ws.ClothoWebSocket;
import org.clothocad.core.util.JSON;
import org.eclipse.jetty.websocket.api.Session;
import org.eclipse.jetty.websocket.api.annotations.OnWebSocketClose;
import org.eclipse.jetty.websocket.api.annotations.OnWebSocketConnect;
import org.eclipse.jetty.websocket.api.annotations.OnWebSocketMessage;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.web.socket.WebSocketHandler;
import org.springframework.web.socket.adapter.jetty.JettyWebSocketHandlerAdapter;
import org.springframework.web.socket.adapter.jetty.JettyWebSocketSession;

/**
 *
 * @author david
 */
public class ClothoJettyHandler extends JettyWebSocketHandlerAdapter {

    private Router router;
    private ClothoWebSocket ws;
    private String id;
    
    private static Logger logger = LoggerFactory.getLogger(ClothoJettyHandler.class);

    public ClothoJettyHandler(WebSocketHandler webSocketHandler, JettyWebSocketSession wsSession) {
        super(webSocketHandler, wsSession);
        id = wsSession.getId();
    }

    public void setRouter(Router router) {
        this.router = router;
        ws = new ClothoWebSocket(id, router);
    }

    @OnWebSocketClose
    public void onClose(int closeCode, String message) {
        ws.onClose(closeCode, message);
    }

    public void send(Message msg) {
        ws.send(msg);

    }

    @OnWebSocketMessage
    public void onMessage(String messageString) {
        try {
            logger.debug("Handing message over to ClothoWebSocket: " + messageString);
            ws.onMessage(messageString);
        }
        catch(Throwable ex)
        {
            logger.debug(this.getClass().getName() + ": onMessageException - " + ex.getMessage());
        }
    }

    public boolean isOpen() {
        return ws.isOpen();
    }

    @OnWebSocketConnect
    public void onOpen(Session session) {
        ws.onOpen(session);

    }

}

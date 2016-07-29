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

/**
 *
 * @author david
 */
@Component
public class ClothoWSHandler implements WebSocketHandler {

    private Router router;

    private ServerSideAPI api;

    private static Logger logger = LoggerFactory.getLogger(ClothoWSHandler.class);

    public void setRouter(Router router) {
        this.router = router;
    }

    @Override
    public void afterConnectionEstablished(WebSocketSession session) throws Exception {
        logger.debug("Session: " + session.getId());
        logger.debug("Session: " + session.getId());
        logger.debug("Session: " + session.getId());
        logger.debug("Session: " + session.getId());

    }

    @Override
    public void handleMessage(WebSocketSession session, WebSocketMessage<?> message) throws Exception {
        logger.debug(message.getPayload().toString());
        logger.debug(message.getPayload().toString());
        logger.debug(message.getPayload().toString());
        logger.debug(message.getPayload().toString());
        logger.debug(message.getPayload().toString());
        logger.debug(message.getPayload().toString());
        logger.debug(message.getPayload().toString());

        logger.debug("" + session.isOpen());

        ClothoWebSocket ws = new ClothoWebSocket(session.getId(), router);

        logger.debug("" + session.isOpen());

        String payload = message.getPayload().toString();

        JSONObject obj = new JSONObject(payload);
        logger.debug("" + session.isOpen());
        String channel = obj.getString("channel");
        Channel realChannel = Channel.alert;

        for (Channel n : Channel.values()) {
            if (n.name().equalsIgnoreCase(channel)) {
                realChannel = n;
            }
        }
        if (realChannel == Channel.alert) {
            realChannel = Channel.autocomplete;
        }

        logger.debug(realChannel.name());
        logger.debug(realChannel.name());
        logger.debug(realChannel.name());
        logger.debug(realChannel.name());
        logger.debug(realChannel.name());
        logger.debug(realChannel.name());
        logger.debug(realChannel.name());
        logger.debug(realChannel.name());

        Object data = obj.get("data");
        String rqID = obj.getString("requestId");
        logger.debug("" + session.isOpen());
//        Map<MessageOption, Object> opt = (Map<MessageOption, Object>) obj.get("options");

        logger.debug("" + session.isOpen());

        Message msg = new Message(realChannel, data, rqID);

        logger.debug(msg.toString());
        logger.debug("" + session.isOpen());
        //Gets trapped and is sad. :<
        router.receiveMessage(ws, msg);
        logger.debug("" + session.isOpen());
        logger.debug("Router is receiving message");

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

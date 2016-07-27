/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.webserver.jetty;

import javax.inject.Inject;
import org.clothocad.core.communication.Router;
import org.clothocad.core.communication.ws.ClothoWebSocket;
import org.eclipse.jetty.websocket.api.annotations.WebSocket;
import org.eclipse.jetty.websocket.servlet.ServletUpgradeRequest;
import org.eclipse.jetty.websocket.servlet.ServletUpgradeResponse;
import org.eclipse.jetty.websocket.servlet.WebSocketCreator;
import org.eclipse.jetty.websocket.servlet.WebSocketServlet;
import org.eclipse.jetty.websocket.servlet.WebSocketServletFactory;

/**
 *
 * @author david
 */
@WebSocket
@SuppressWarnings("serial")
public class ClothoServlet extends WebSocketServlet {

    public Router router;
    
    public ClothoServlet(final Router router) {
        this.router = router;//To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void configure(WebSocketServletFactory factory) {
        factory.getPolicy().setIdleTimeout(10000);
        factory.setCreator(new WebSocketCreator() {
            @Override
            public Object createWebSocket(ServletUpgradeRequest sur, ServletUpgradeResponse sur1) {

                return new ClothoWebSocket(sur.getHttpServletRequest().getSession().getId(), router);
            }
        });
    }
}

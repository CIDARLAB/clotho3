package org.clothocad.webserver.jetty;

import javax.inject.Named;
import javax.servlet.http.HttpServletRequest;

import org.clothocad.core.layers.communication.connection.ws.ClothoWebSocket;
import org.eclipse.jetty.server.Handler;
import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.server.handler.ContextHandler;
import org.eclipse.jetty.server.handler.HandlerList;
import org.eclipse.jetty.server.handler.ResourceHandler;
import org.eclipse.jetty.websocket.WebSocket;
import org.eclipse.jetty.websocket.WebSocketHandler;
import org.eclipse.jetty.server.session.SessionHandler;

public class ClothoWebserver {

	public ClothoWebserver(@Named("port") int nPort) 
			throws Exception {
		
		Server server = new Server(nPort);
 
        /** for WEB resources (HTML, CSS, JavaScript etc.) **/
        ResourceHandler resHandler = new ResourceHandler();
        resHandler.setDirectoriesListed(true);
        
        //Regardless of repository, the URI to access is http://localhost:8080/#/
        //Set directory to clotho3-web (this repository)
        resHandler.setWelcomeFiles(new String[]{ "index.html" }); 
        resHandler.setResourceBase("./clotho3-web/");
        
        //Set directory to Max's Github repository
//        resHandler.setWelcomeFiles(new String[]{ "index.html" }); 
//        resHandler.setResourceBase("./../clotho_multiMode/app/");

        /** Clotho3.0 Java Websocket **/
        System.out.println("Ernst (not for sb6.0 demo, but eventually), this should be pulling from DB namespaced by the View's uuid to avoid collisions in naming.  Not in flatfiles");
        WebSocketHandler wsHandler = new WebSocketHandler() {
        	@Override
        	public WebSocket doWebSocketConnect(HttpServletRequest request, String protocol) {
        		return new ClothoWebSocket(request.getSession().getId());
        	}
        };        
        ContextHandler contextHandler = new ContextHandler();
        contextHandler.setContextPath("/websocket");
        contextHandler.setHandler(wsHandler);
        
        /** Session Handling ***/
        SessionHandler sessionHandler = new SessionHandler();
        
        HandlerList handlers = new HandlerList();
        handlers.setHandlers(new Handler[] { 
        		sessionHandler, 
        		resHandler,        		
        		wsHandler });
        server.setHandler(handlers);
 
        server.start();
        server.join();
        
        //System.out.println("Ernst, please silence these println statements");
        
        //System.out.println("The Clotho Webserver is running...");
    }
	
/**
	public static void main(String[] args) {
		try {
			new ClothoWebserver();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
 **/
}
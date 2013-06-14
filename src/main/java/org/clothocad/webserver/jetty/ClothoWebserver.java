package org.clothocad.webserver.jetty;

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

	public ClothoWebserver() 
			throws Exception {
		
		Server server = new Server(8080);
 
        /** for WEB resources (HTML, CSS, JavaScript etc.) **/
        System.out.println("Ernst, this needs to be changed to fetch things from the database instead of flatfiles.  The organization should use the UUID of the View to namespace things and avoid naming collisions  between projects.");
        ResourceHandler resHandler = new ResourceHandler();
        resHandler.setDirectoriesListed(true);
        resHandler.setWelcomeFiles(new String[]{ "index.html" }); 
        resHandler.setResourceBase("./clotho3-web/");
        
        /** Clotho3.0 Java Websocket **/
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
        
        System.out.println("Ernst, Please silence the verbose logging of Jetty");        
        
        System.out.println("The Clotho Webserver is running...");
    }
	
	public static void main(String[] args) {
		try {
			new ClothoWebserver();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
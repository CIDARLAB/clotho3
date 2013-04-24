package org.clothocad.webserver.jetty.rest;

import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.eclipse.jetty.util.log.Log;
import org.eclipse.jetty.websocket.WebSocket;
import org.eclipse.jetty.websocket.WebSocketServlet;

public class ClothoRestServlet 
		extends WebSocketServlet {

	private final Set<ClothoWebSocket> _members = 
			new HashSet<ClothoWebSocket>();
	
	@Override
	protected void doGet(
			HttpServletRequest request, 
			HttpServletResponse response)
			throws ServletException, IOException {
		getServletContext().getNamedDispatcher("default").forward(request,response);    
	}
	
	@Override
	public WebSocket doWebSocketConnect(
			HttpServletRequest request, 
			String protocol) {
		return new ClothoWebSocket();
	}
	
}

package org.clothocad.core.jetty;

import java.io.IOException;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.eclipse.jetty.server.Request;
import org.eclipse.jetty.server.handler.AbstractHandler;

public class RESTHandler 
	extends AbstractHandler {
	
	public void handle(String target, Request baseRequest, HttpServletRequest request, HttpServletResponse response)
	    throws IOException, ServletException {
		
	    String username = baseRequest.getParameter("username");
	    String password = baseRequest.getParameter("password");

	    // TODO: forward the incoming request to the Clotho Core
	    
	    // write the response
	    response.setContentType("text/html;charset=utf-8");
	    response.setStatus(HttpServletResponse.SC_OK);
	    baseRequest.setHandled(true);	
	    response.getWriter().println("Hello ..."+username);
	}
}

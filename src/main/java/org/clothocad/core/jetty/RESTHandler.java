package org.clothocad.core.jetty;

import java.io.IOException;
import java.util.Enumeration;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.eclipse.jetty.server.Request;
import org.eclipse.jetty.server.handler.AbstractHandler;

public class RESTHandler 
	extends AbstractHandler {
	
	public void handle(String target, Request baseRequest, HttpServletRequest request, HttpServletResponse response)
	    throws IOException, ServletException {
		
		// now, let's print the request's parameters
		Enumeration<String> paramNames = baseRequest.getParameterNames();
		while(paramNames.hasMoreElements()) {
			String sParamName = paramNames.nextElement();
			System.out.println("[RESTHandler.handle] -> "+sParamName+" -> "+baseRequest.getParameter(sParamName));
		}
		
	    // TODO: forward the incoming request to the Clotho Core
	    
	    // write the response
	    response.setContentType("text/html;charset=utf-8");
	    response.setStatus(HttpServletResponse.SC_OK);
	    baseRequest.setHandled(true);	
	    response.getWriter().println("This is the server's response");
	}
}

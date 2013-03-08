package org.clothocad.core.jetty;

import java.io.IOException;
import java.util.Enumeration;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.eclipse.jetty.server.Request;
import org.eclipse.jetty.server.handler.AbstractHandler;

public class ClothoDemoHandler 
		extends AbstractHandler {
	
	public void handle(String target, Request baseRequest, HttpServletRequest request, HttpServletResponse response) 
	        throws IOException, ServletException {

		System.out.println("****** REQUEST INFORMATION ****");
		// print the information
		System.out.println("target: "+target);

		System.out.println("*** BASE REQUEST ***");
		
		System.out.println("*** PARAMETERS");
		Enumeration<String> paramNames = baseRequest.getParameterNames();
		while(paramNames.hasMoreElements()) {
			String paramName = paramNames.nextElement();
			System.out.println(paramName+": "+baseRequest.getParameter(paramName));			
		} 
		
		System.out.println("*** SERVLET REQUEST ***");
		
		
		response.setContentType("text/html;charset=utf-8");
		response.setStatus(HttpServletResponse.SC_OK);
		
		baseRequest.setHandled(true);
		
		response.getWriter().println("<h1>Hello World</h1>");
	}
	
}

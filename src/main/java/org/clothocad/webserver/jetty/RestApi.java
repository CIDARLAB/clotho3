package org.clothocad.webserver.jetty;

import java.io.IOException;
import java.util.Map;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

@SuppressWarnings("serial")
public class RestApi extends HttpServlet {
    protected void doGet(HttpServletRequest request, 
    	HttpServletResponse response) 
    throws ServletException, IOException {

    	response.setContentType("text/html");

    	String path = request.getPathInfo();
    	Map<String, String[]> params = request.getParameterMap();

		response.setStatus(HttpServletResponse.SC_OK);
		response.getWriter().println("<h1>Buongiorno Mijo</h1>");
		response.getWriter().println("path: " + path + "<br>");
		response.getWriter().println("params: " + params);
    }
}

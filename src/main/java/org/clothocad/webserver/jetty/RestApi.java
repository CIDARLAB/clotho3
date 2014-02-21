package org.clothocad.webserver.jetty;

import java.io.IOException;
import java.util.Map;
import java.util.Enumeration;

import org.clothocad.core.communication.*;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

@SuppressWarnings("serial")
public class RestApi extends HttpServlet {

    protected void doGet(HttpServletRequest request, 
    	HttpServletResponse response) 
    throws ServletException, IOException {

    	response.setContentType("text/json");

    	String path = request.getPathInfo();
        String paramValue;

        switch (path) {
            // Need to add possibility of having trailing slash.
            case "/autocomplete/":
            case "/autocomplete":
                paramValue = (String) request.getParameter("userText");
                if (paramValue != null) {
                    response.setStatus(HttpServletResponse.SC_OK);
                    response.getWriter().println("<h1>autocomplete</h1>");
                } else {
                    response.sendError(response.SC_BAD_REQUEST, "Parameter userText required");
                }
                break;
            case "/autocompleteDetail/":
            case "/autocompleteDetail":
                paramValue = (String) request.getParameter("uuid");
                if (paramValue != null) {
                    response.setStatus(HttpServletResponse.SC_OK);
                    response.getWriter().println("<h1>autocompleteDetail</h1>");
                } else {
                    response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
                }
                break;
        }
    }
}

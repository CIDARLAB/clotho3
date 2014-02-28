package org.clothocad.webserver.jetty;

import java.io.IOException;
import java.util.Map;
import java.util.Enumeration;

import com.google.inject.Guice;
import com.google.inject.Injector;

import org.clothocad.core.communication.*;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.execution.Mind;
import org.clothocad.core.persistence.mongodb.MongoDBModule;
// import org.clothocad.core.testers.ClothoTestModule;
// This is in the test package

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

@SuppressWarnings("serial")
public class RestApi extends HttpServlet {

    private static ServerSideAPI api;
    private static ServerSideAPI unprivilegedUser;
    private static Persistor persistor;
    private static Mind mind;
    private static Router router;

    private void createAPI() {
        // Injector injector = Guice.createInjector(new ClothoTestModule(), new MongoDBModule());
        // persistor = injector.getInstance(Persistor.class);
        // router = injector.getInstance(Router.class);
        // mind = new Mind();
        // api = new ServerSideAPI(mind, persistor, router, null);
        // persistor.connect();
        // mind.setConnection(new TestConnection("test"));
    }

    protected void doGet(HttpServletRequest request, 
    	HttpServletResponse response) throws ServletException, IOException {

    	response.setContentType("text/json");

    	String path = request.getPathInfo();

        switch (path) {
            // Needed to add possibility of having trailing slash.
            case "/autocomplete/":
            case "/autocomplete":
                String paramUserText = (String) request.getParameter("userText");
                if (paramUserText != null) {
                    response.setStatus(HttpServletResponse.SC_OK);
                    response.getWriter().println("<h1>autocomplete</h1>");
                } else {
                    response.sendError(response.SC_BAD_REQUEST, "Parameter userText required");
                }
                break;
            case "/autocompleteDetail/":
            case "/autocompleteDetail":
                String paramUUID = (String) request.getParameter("uuid");
                if (paramUUID != null) {
                    response.setStatus(HttpServletResponse.SC_OK);
                    response.getWriter().println("<h1>autocompleteDetail</h1>");
                } else {
                    response.sendError(response.SC_BAD_REQUEST, "Parameter uuid required");
                }
                break;
            case "/submit/":
            case "/submit":
                String paramCommand = (String) request.getParameter("command");
                if (paramCommand != null) {
                    response.setStatus(HttpServletResponse.SC_OK);
                    response.getWriter().println("<h1>submit</h1>");
                } else {
                    response.sendError(response.SC_BAD_REQUEST, "Parameter command required");
                }
                break;
            case "/learn/":
            case "/learn":
                // Need to deal with other case of learn where an object is passed in to SSA
                String paramUI = (String) request.getParameter("userInput");
                String paramCommand1 = (String) request.getParameter("command");
                if (paramUI != null && paramCommand1 != null) {
                    response.setStatus(HttpServletResponse.SC_OK);
                    response.getWriter().println("<h1>learn</h1>");
                } else {
                    response.sendError(response.SC_BAD_REQUEST, "Parameter command & userInput required");
                }
                break;
            case "/login/":
            case "/login":
                String paramUName = (String) request.getParameter("username");
                String paramPass = (String) request.getParameter("password");
                if (paramUName != null && paramPass != null) {
                    response.setStatus(HttpServletResponse.SC_OK);
                    response.getWriter().println("<h1>login</h1>");
                } else {
                    response.sendError(response.SC_BAD_REQUEST, "Parameter username & password required");
                }
                break;
        }
    }
}

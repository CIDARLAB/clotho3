package org.clothocad.webserver.jetty;

import java.io.IOException;
import java.util.Map;

import com.google.inject.Guice;
import com.google.inject.Injector;

import org.clothocad.core.communication.*;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.execution.Mind;
import org.clothocad.core.persistence.mongodb.MongoDBModule;
import org.apache.shiro.SecurityUtils;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

@SuppressWarnings("serial")
public class RestApi extends HttpServlet {

    private static ServerSideAPI api;
    private static Persistor persistor;
    private static Mind mind;
    private static Router router;

    // retrieve mind object/create new one - session id w/ shiro
    // minds are in router

    // http://shiro.apache.org/
    // 2 methods:
    // pass request to serversideapi w/ mind obj
    // get : shareable through uuid
    // post : change value of shareable with new val

    protected void doGet(HttpServletRequest request, 
    	HttpServletResponse response) throws ServletException, IOException {

    	response.setContentType("text/json");

    	String path = request.getPathInfo();
        String id = path.split("/")[0];

        if (id == null) {
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
            response.getWriter().println("requestId is required - the proper GET url format is .../rest/{id}");
        } else if (SecurityUtils.getSubject().isAuthenticated()) { // the user is logged in. Will this ever evaluate to true? Not sure how this works exactly.
            String username = SecurityUtils.getSubject().getPrincipal().toString();
            mind.setUsername(username);
            persistor.save(mind);
            // How do we retrieve/generate router?
            // What is the requestId? Does this work?
            api = new ServerSideAPI(mind, persistor, router, id);
            // What is the object we will pass into api.get(Object o) here? 
            // How can it be represented (i.e. is there an id we can use to retrieve it? How would it work otherwise over http?)
            Map<String, Object> result = api.get();

            // here I want to generate a json response mapping string -> object

            response.setStatus(HttpServletResponse.SC_OK);
            response.getWriter().println("You are logged in");

        } else {
            // What modules would you use in createInjector if we are using unlogged in mind? 
            // Or will we disallow this?
            // Injector injector = Guice.createInjector(new ClothoTestModule(), new MongoDBModule());
            persistor = injector.getInstance(Persistor.class);
            router = injector.getInstance(Router.class);
            mind = new Mind();
            api = new ServerSideAPI(mind, persistor, router, null);
            persistor.connect();
            // mind.setConnection(new TestConnection("test"));

            response.setStatus(HttpServletResponse.SC_OK);
            response.getWriter().println("You are not logged in");
        }
    }

    protected void doPost(HttpServletRequest request, 
        HttpServletResponse response) throws ServletException, IOException {

        response.setContentType("text/json");

        String path = request.getPathInfo();
        String id = path.split("/")[0];

        if (id == null) {
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
            response.getWriter().println("requestId is required - the proper POST url format is .../rest/{id}");
        } else if (SecurityUtils.getSubject().isAuthenticated()) {
            String username = SecurityUtils.getSubject().getPrincipal().toString();
            mind.setUsername(username);
            persistor.save(mind);
            api = new ServerSideAPI(mind, persistor, router, id);
            // again, how do we pass an object in http? Does it have an id I can grab?
            ObjectId response = api.set();
            response.setStatus(HttpServletResponse.SC_OK);
            response.getWriter().println("Logged in to post");
        }
    }
}

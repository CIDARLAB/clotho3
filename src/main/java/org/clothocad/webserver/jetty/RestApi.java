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
// import org.clothocad.core.testers.ClothoTestModule;
// This is in the test package

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

        if (path == "get/trial/") {
            if (SecurityUtils.getSubject().isAuthenticated()) {
                String username = SecurityUtils.getSubject().getPrincipal().toString();
                mind.setUsername(username);
                persistor.save(mind);
                response.setStatus(HttpServletResponse.SC_OK);
                response.getWriter().println("You are logged in");
            } else {
                response.setStatus(HttpServletResponse.SC_OK);
                response.getWriter().println("You are not logged in");
            }
        }
    }

    protected void doPost(HttpServletRequest request, 
        HttpServletResponse response) throws ServletException, IOException {

        response.setContentType("text/json");

        String path = request.getPathInfo();

        if (path == "post/trial/") {
            if (SecurityUtils.getSubject().isAuthenticated()) {
                String username = SecurityUtils.getSubject().getPrincipal().toString();
                mind.setUsername(username);
                persistor.save(mind);
                response.setStatus(HttpServletResponse.SC_OK);
                response.getWriter().println("<h1>autocomplete</h1>");
            }
        }
    }
}

package org.clothocad.webserver.jetty;

import java.io.IOException;
import java.util.Map;
import java.util.HashMap;

import com.google.inject.Guice;
import com.google.inject.Injector;
import javax.inject.Named;
import javax.inject.Inject;

import org.clothocad.core.communication.*;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.execution.Mind;
import org.clothocad.core.persistence.mongodb.MongoDBModule;
import org.clothocad.core.util.JSON;
import org.clothocad.core.datums.ObjBase;

import org.apache.shiro.subject.Subject;
import org.apache.shiro.SecurityUtils;
import org.bson.types.ObjectId;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

@SuppressWarnings("serial")
public class RestApi extends HttpServlet {

    private static Router router;
    private static Message m;
    private static RestConnection rc = new RestConnection("RestConnection");;
    private static Map<String, Object> messageMap;

    // retrieve mind object/create new one - session id w/ shiro
    // minds are in router

    // http://shiro.apache.org/
    // 2 methods:
    // pass request to serversideapi w/ mind obj
    // get : shareable through uuid
    // post : change value of shareable with new val

    public RestApi(Router router) {
        this.router = router;
    }


    protected void doGet(HttpServletRequest request, 
    	HttpServletResponse response) throws ServletException, IOException {

    	response.setContentType("application/json");

    	String id = request.getPathInfo().split("/")[1];

        // We build our new message
        messageMap = new HashMap<String, Object>();
        messageMap.put("data", id);
        messageMap.put("channel", "get");

        m = new Message(messageMap);

        // Now we send that message to the router
        this.router.receiveMessage(this.rc, m);

        String result = this.rc.getResult().toString();
        response.getWriter().write(result);

        if (!result.contains("FAILURE")) {
            response.setStatus(HttpServletResponse.SC_OK);
        } else {
            response.setStatus(HttpServletResponse.SC_NOT_FOUND);
        }
    }

    protected void doPost(HttpServletRequest request, 
        HttpServletResponse response) throws ServletException, IOException {

        response.setContentType("application/json");

        // Need to add checks for id
        String id = request.getPathInfo().split("/")[1];
        // For set, we need other parameters. What are they?

        // We build our new message
        messageMap = new HashMap<String, Object>();
        messageMap.put("data", id);
        messageMap.put("channel", "set");

        m = new Message(messageMap);

        // Now we send that message to the router
        this.router.receiveMessage(this.rc, m);

        // Get the result & check to see if it was successful/if it failed
        String result = this.rc.getResult().toString();
        response.getWriter().write(result);

        if (!result.contains("FAILURE")) {
            response.setStatus(HttpServletResponse.SC_OK);
        } else {
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
        }
    }
}

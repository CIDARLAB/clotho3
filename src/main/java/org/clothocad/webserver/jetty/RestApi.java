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

    // retrieve mind object/create new one - session id w/ shiro
    // minds are in router

    // http://shiro.apache.org/
    // 2 methods:
    // pass request to serversideapi w/ mind obj
    // get : shareable through uuid
    // post : change value of shareable with new val

    @Inject
    public void RestApi(Router router) {
        this.router = router;
    }


    protected void doGet(HttpServletRequest request, 
    	HttpServletResponse response) throws ServletException, IOException {

    	response.setContentType("text/json");

    	String path = request.getPathInfo();
        String id = path.split("/")[1];

        Subject subject = SecurityUtils.getSubject();

        if (id.equals("")) { // no id has been supplied
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
            response.getWriter().println("requestId is required - the proper GET url format is .../rest/{id}");
        } else {
            // We build our new message
            Map<String, Object> messageMap = new HashMap<String, Object>();
            messageMap.put("data", id);
            messageMap.put("channel", "get");
    
            m = new Message(messageMap);

            // Now we send that message to the router
            System.out.println(this.router);
            this.router.receiveMessage(this.rc, m);

            response.getWriter().write(this.rc.getResult().toString());
            response.setStatus(HttpServletResponse.SC_OK);
        }
    }

    protected void doPost(HttpServletRequest request, 
        HttpServletResponse response) throws ServletException, IOException {

        response.setContentType("text/json");

        String path = request.getPathInfo();
        String id = path.split("/")[0];

        if (id.equals("")) { // no id has been supplied
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
            response.getWriter().println("requestId is required - the proper POST url format is .../rest/{id}");
        } else if (SecurityUtils.getSubject().isAuthenticated()) {
            // String username = SecurityUtils.getSubject().getPrincipal().toString();
            // mind.setUsername(username);
            // persistor.save(mind);
            // api = new ServerSideAPI(mind, persistor, null, null);
            // again, how do we pass an object in http? Does it have an id I can grab?
            // ObjectId response = api.set();
            response.setStatus(HttpServletResponse.SC_OK);
            response.getWriter().println("Logged in to post");
        }
    }
}

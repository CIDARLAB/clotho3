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

    private static ServerSideAPI api;
    private static Persistor persistor;
    private static Mind mind;
    private static Router router;
    private static Message m;

    // retrieve mind object/create new one - session id w/ shiro
    // minds are in router

    // http://shiro.apache.org/
    // 2 methods:
    // pass request to serversideapi w/ mind obj
    // get : shareable through uuid
    // post : change value of shareable with new val

    @Inject
    public void RestApi(@Named("apiRouter") Router apiRouter,
        @Named("persistorObj") Persistor persistorObj) {
        this.persistor = persistorObj;
        this.router = apiRouter;
    }


    protected void doGet(HttpServletRequest request, 
    	HttpServletResponse response) throws ServletException, IOException {
        // response.getWriter().write(json.toString());

    	response.setContentType("text/json");

    	String path = request.getPathInfo();
        String id = path.split("/")[1];

        Subject subject = SecurityUtils.getSubject();

        if (id.equals("")) { // no id has been supplied
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
            response.getWriter().println("requestId is required - the proper GET url format is .../rest/{id}");
        } else if (subject.isAuthenticated()) { // the user is logged in.
            // We can do this all through router api - apiRouter.receiveMessage(...)

            mind = getAuthenticatedMind(subject.getPrincipal().toString());
            // how do we run mind.setconnection() ?
            persistor.save(mind);
            // How do we retrieve/generate router?
            api = new ServerSideAPI(mind, this.persistor, this.router, null);

            // In order for this to work properly, id has to be incapsulated in a Message object
            // Map<String, Object> result = api.receiveMessage(id);
            // String jsonResult = JSON.serialize(result);

            response.setStatus(HttpServletResponse.SC_OK);
            // response.getWriter().println(jsonResult);
        } else {
            // router.receiveMessage(..) requires connection & Message request
            // Map<String, Object> result = apiRouter.receiveMessage(id);
            // String jsonResult = JSON.serialize(result);
            mind = new Mind();
            // how do we run mind.setconnection() ?
            // persistor.save(mind);

            api = new ServerSideAPI(mind, this.persistor, this.router, null);

            Map<String, Object> result = api.get(id);
            System.out.println("hey");
            String jsonResult = JSON.serialize(result);

            response.setStatus(HttpServletResponse.SC_OK);
            response.getWriter().write(jsonResult.toString());
        }


        // {
        //     Map<String, Object> messageMap = new HashMap<String, Object>();
        //     messageMap.put("requestId", null);
        //     messageMap.put("data", id);
        //     m = new Message(messageMap);
        //     apiRouter.receiveMessage(null, m);
        // }
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
            String username = SecurityUtils.getSubject().getPrincipal().toString();
            mind.setUsername(username);
            persistor.save(mind);
            api = new ServerSideAPI(mind, persistor, null, null);
            // again, how do we pass an object in http? Does it have an id I can grab?
            // ObjectId response = api.set();
            response.setStatus(HttpServletResponse.SC_OK);
            response.getWriter().println("Logged in to post");
        }
    }


    private Mind getAuthenticatedMind(String username)  {
        //XXX: this whole method is janky
        Map<String,Object> query = new HashMap();
        query.put("username", username);
        query.put("className", Mind.class.getCanonicalName());
        try {
            Iterable<ObjBase> minds = persistor.find(query);
            Mind mind;
            
            if (!minds.iterator().hasNext()){
                mind = new Mind();
            } else {
                mind = (Mind) minds.iterator().next();
            }

            return mind;
        } catch (Exception ex){
            ex.printStackTrace();
            throw ex;
        }
    }
}

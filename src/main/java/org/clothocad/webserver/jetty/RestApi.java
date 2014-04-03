package org.clothocad.webserver.jetty;

import java.io.IOException;
import java.io.BufferedReader;
import java.util.Map;
import java.util.HashMap;
import java.util.Collection;

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
import javax.servlet.http.Part;

/**
 *
 * @author mosnicholas
 */
@SuppressWarnings("serial")
public class RestApi extends HttpServlet {

    private static Router router;
    private static Message m;
    private static RestConnection rc = new RestConnection("RestConnection");
    // TODO:
    // Set
    // Test set with following url , change sequence: https://localhost:8443/rest/52410adf50763ce31f941915

    // http://aaronparecki.com/articles/2012/07/29/1/oauth2-simplified
    // http://stackoverflow.com/questions/15051712/how-to-do-authentication-with-a-rest-api-right-browser-native-clients

    // write create api endpoint
    // request body json for set/create
    // look into http authentication

    // http://shiro.apache.org/
    // 2 methods:
    // pass request to serversideapi w/ mind obj
    // post : change value of shareable with new val
    public RestApi(Router router) {
        this.router = router;
    }


    protected void doGet(HttpServletRequest request, 
    	HttpServletResponse response) throws ServletException, IOException {

    	response.setContentType("application/json");

    	String id = request.getPathInfo().split("/")[1];
        
        // We build our new message
        m = new Message(Channel.get, id, null, null);

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

    protected void doDelete(HttpServletRequest request, 
        HttpServletResponse response) throws ServletException, IOException {

        response.setContentType("application/json");

        String id = request.getPathInfo().split("/")[1];

        // We build our new message
        m = new Message(Channel.destroy, id, null, null);

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
        // Params for set: Id of object, new fields you want to update
        // Put new parameters as json in the request body.

        response.setContentType("application/json");

        // Need to add checks for id
        String id = request.getPathInfo().split("/")[1];
        
        Map<String, String> p = getRequestBody(request.getReader());
        p.put("id", id);

        // We build our new message
        m = new Message(Channel.set, p, null, null);

        // Now we send that message to the router
        this.router.receiveMessage(this.rc, m);

        // Get the result & check to see if it was successful/if it failed
        String result = this.rc.getResult().toString();;
        
        response.getWriter().write(result);

        if (!result.contains("FAILURE")) {
            response.setStatus(HttpServletResponse.SC_OK);
        } else {
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
        }
    }

    private Map<String, String> getRequestBody(BufferedReader reader) {
        String parts[];
        String keyValue[] = new String[2];
        String currChar = "";
        String key = "";

        Map<String, String> map = new HashMap<String, String>();

        try {
            parts = reader.readLine().split("&");
            for (String kv : parts) {
                keyValue = kv.split("=");
                map.put(keyValue[0], keyValue[1]);
            }
        } catch (IOException ie) {
            map = new HashMap<String, String>();
            map.put("IOException", ie.toString());
        }

        return map;
    }
}

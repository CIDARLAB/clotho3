package org.clothocad.webserver.jetty;

import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.io.BufferedReader;
import java.util.Map;
import java.util.HashMap;

import org.clothocad.core.communication.*;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.execution.Mind;
import org.clothocad.core.util.JSON;

import org.apache.commons.codec.binary.Base64;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

/**
 *
 * @author mosnicholas
 */
@SuppressWarnings("serial")
public class RestApi extends HttpServlet {

    private static Router router;
    private static Message m;
    private static Message loginMessage, logoutMessage;
    private static Map<String, String> loginMap;
    private static RestConnection rc = new RestConnection("RestConnection");
    // Test set with following url , change sequence: https://localhost:8443/rest/52410adf50763ce31f941915

    // http://stackoverflow.com/questions/15051712/how-to-do-authentication-with-a-rest-api-right-browser-native-clients
    // http://shiro-user.582556.n2.nabble.com/Shiro-and-RESTful-web-services-td5539212.html
    // http://stackoverflow.com/questions/319530/restful-authentication?rq=1
    // http://stackoverflow.com/questions/454355/security-of-rest-authentication-schemes

    // should i login through router & then logout?

    // look into http authentication

    // http://shiro.apache.org/
    public RestApi(Router router) {
        this.router = router;
    }

    protected void doGet(HttpServletRequest request, 
    	HttpServletResponse response) throws ServletException, IOException {

    	response.setContentType("application/json");

        String[] unamePass = getBasicAuth(request.getHeader("Authorization"));

        if (unamePass != null) {
            loginMap = new HashMap<String, String>();
            loginMap.put("username", unamePass[0]);
            loginMap.put("password", unamePass[1]);
            loginMessage = new Message(Channel.login, loginMap, null, null);
            this.router.receiveMessage(this.rc, loginMessage);
        }

    	String id = request.getPathInfo().split("/")[1];

        // We build our new message
        m = new Message(Channel.get, id, null, null);

        // Now we send that message to the router
        this.router.receiveMessage(this.rc, m);

        String result = this.rc.getResult().toString();

        if (unamePass != null) {
            logoutMessage = new Message(Channel.logout, loginMap, null, null);
            this.router.receiveMessage(this.rc, logoutMessage);
        }
        
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

        String[] unamePass = getBasicAuth(request.getHeader("Authorization"));

        if (unamePass != null) {
            loginMap = new HashMap<String, String>();
            loginMap.put("username", unamePass[0]);
            loginMap.put("password", unamePass[1]);
            loginMessage = new Message(Channel.login, loginMap, null, null);
            this.router.receiveMessage(this.rc, loginMessage);
        }

        String id = request.getPathInfo().split("/")[1];

        // We build our new message
        m = new Message(Channel.destroy, id, null, null);

        // Now we send that message to the router
        this.router.receiveMessage(this.rc, m);

        String result = this.rc.getResult().toString();

        if (unamePass != null) {
            logoutMessage = new Message(Channel.logout, loginMap, null, null);
            this.router.receiveMessage(this.rc, logoutMessage);
        }

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

        String[] unamePass = getBasicAuth(request.getHeader("Authorization"));

        if (unamePass != null) {
            loginMap = new HashMap<String, String>();
            loginMap.put("username", unamePass[0]);
            loginMap.put("password", unamePass[1]);
            loginMessage = new Message(Channel.login, loginMap, null, null);
            this.router.receiveMessage(this.rc, loginMessage);
        }

        String id = request.getPathInfo().split("/")[1];

        Map<String, String> p = getRequestBody(request.getReader());
        p.put("id", id);

        // We build our new message
        m = new Message(Channel.set, p, null, null);

        // Now we send that message to the router
        this.router.receiveMessage(this.rc, m);

        // Get the result & check to see if it was successful/if it failed
        String result = this.rc.getResult().toString();

        if (unamePass != null) {
            logoutMessage = new Message(Channel.logout, loginMap, null, null);
            this.router.receiveMessage(this.rc, logoutMessage);
        }
        
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

    private String[] getBasicAuth(String authHeader) {
        try {
            String credentials = new String(Base64.decodeBase64(authHeader), "UTF-8");
            return credentials.split(":");
        } catch (UnsupportedEncodingException uee) {
            return null;
        } catch (NullPointerException nee) {
            return null;
        }
    }
}

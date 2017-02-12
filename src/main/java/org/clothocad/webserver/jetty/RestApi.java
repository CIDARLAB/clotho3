package org.clothocad.webserver.jetty;

import org.clothocad.core.communication.Channel;
import org.clothocad.core.communication.Message;
import org.clothocad.core.communication.RestConnection;
import org.clothocad.core.communication.Router;
import org.clothocad.core.persistence.Persistor;


import java.io.BufferedReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import org.apache.shiro.authz.UnauthorizedException;

@SuppressWarnings("serial")
public class RestApi extends HttpServlet {

    /*
    
    REST API servlet is at idontremember.url/data/*
    
    doGet: get, getAll, and query
    doPost: pterry much any request you want, just specify the channel.
    
    haven't modified doDelete or doPut yet.
    
    
     */
    private static Router router;
    private static Persistor persistor;
    private static Message m, loginMessage, logoutMessage;
    private static Map<String, String> loginMap;
    private static RestConnection rc = new RestConnection("RestConnection");

    // http://stackoverflow.com/questions/15051712/how-to-do-authentication-with-a-rest-api-right-browser-native-clients
    // http://shiro-user.582556.n2.nabble.com/Shiro-and-RESTful-web-services-td5539212.html
    // http://stackoverflow.com/questions/319530/restful-authentication?rq=1
    // http://stackoverflow.com/questions/454355/security-of-rest-authentication-schemes
    // http://shiro.apache.org/
    public RestApi(Router router) {
        this.router = router;
    }

    protected void doGet(HttpServletRequest request,
            HttpServletResponse response) throws ServletException, IOException {

        response.addHeader("Access-Control-Allow-Origin", "*");
        response.setContentType("application/json");

        String[] pathID = request.getPathInfo().split("/");
        Map<String, String> body = getRequestBody(request.getReader());

        String[] auth = new String(request.getHeader("Authorization")).split(":");
        
        System.out.println("\n\n\n" + Arrays.toString(auth) + "\n\n\n");
        
        login(auth);
        
        
        String id = pathID[2];
        String data = pathID[3];
        
        
        switch (id) {
                case "getAll":
            }
        
        
        
        
        
        
        try {
            this.router.receiveMessage(this.rc, m);
        } catch (UnauthorizedException ue) {
            response.setStatus(HttpServletResponse.SC_UNAUTHORIZED);
            response.addHeader("WWW-Authenticate", "Basic realm=\"Clotho Rest\"");
            response.addHeader("HTTP/1.0 401", "Unauthorized");
            response.getWriter().write("{\"error\": \"unauthorized access of page\"}");
            return;
        }

        String result = this.rc.getResult().toString();
        System.out.println("\n\n\n" + result + "\n\n\n");
        //logout(unamePass);

        response.getWriter().write(result);

        if (result.contains("FAILURE")) {
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
        } else {
            response.setStatus(HttpServletResponse.SC_OK);
        }
        
        logout(auth);
        
//        System.out.println("\n\n\n" + body + "\n\n\n");

    }

    protected void doDelete(HttpServletRequest request,
            HttpServletResponse response) throws ServletException, IOException {
        response.addHeader("Access-Control-Allow-Origin", "*");
        response.setContentType("application/json");

        String[] pathID = request.getPathInfo().split("/");
        Map<String, String> body = getRequestBody(request.getReader());
        
        String[] auth = new String(request.getHeader("Authorization")).split(":");
        
        login(auth);
        
        try {
            this.router.receiveMessage(this.rc, m);
        } catch (UnauthorizedException ue) {
            response.setStatus(HttpServletResponse.SC_UNAUTHORIZED);
            response.addHeader("WWW-Authenticate", "Basic realm=\"Clotho Rest\"");
            response.addHeader("HTTP/1.0 401", "Unauthorized");
            response.getWriter().write("{\"error\": \"unauthorized access of page\"}");
            return;
        }

        String result = this.rc.getResult().toString();
        System.out.println("\n\n\n" + result + "\n\n\n");
        //logout(unamePass);

        response.getWriter().write(result);

        if (result.contains("FAILURE")) {
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
        } else {
            response.setStatus(HttpServletResponse.SC_OK);
        }
        
        logout(auth);
    }

    protected void doPost(HttpServletRequest request,
            HttpServletResponse response) throws ServletException, IOException {
        response.addHeader("Access-Control-Allow-Origin", "*");
        response.setContentType("application/json");

        String[] pathID = request.getPathInfo().split("/");
        Map<String, String> body = getRequestBody(request.getReader());
        
        String[] auth = new String(request.getHeader("Authorization")).split(":");
        
        String id = pathID[2];
        
        if (id.equals("createUser")) {
            Map<String,String> credentials = new HashMap<>();
            credentials.put("username", auth[0]);
            credentials.put("credentials", auth[1]);
            credentials.put("displayname", auth[0]);
            m = new Message(Channel.createUser, credentials, null, null);
        }

        login(auth);
        
        try {
            this.router.receiveMessage(this.rc, m);
        } catch (UnauthorizedException ue) {
            response.setStatus(HttpServletResponse.SC_UNAUTHORIZED);
            response.addHeader("WWW-Authenticate", "Basic realm=\"Clotho Rest\"");
            response.addHeader("HTTP/1.0 401", "Unauthorized");
            response.getWriter().write("{\"error\": \"unauthorized access of page\"}");
            return;
        }

        String result = this.rc.getResult().toString();
        System.out.println("\n\n\n" + result + "\n\n\n");
        //logout(unamePass);

        response.getWriter().write(result);

        if (result.contains("FAILURE")) {
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
        } else {
            response.setStatus(HttpServletResponse.SC_OK);
        }
        
        logout(auth);
    }
        
                        

    protected void doPut(HttpServletRequest request,
            HttpServletResponse response) throws ServletException, IOException {
        response.addHeader("Access-Control-Allow-Origin", "*");
        response.setContentType("application/json");

        String[] pathID = request.getPathInfo().split("/");
        Map<String, String> body = getRequestBody(request.getReader());
        
        String[] auth = new String(request.getHeader("Authorization")).split(":");
        
        login(auth);
        
        try {
            this.router.receiveMessage(this.rc, m);
        } catch (UnauthorizedException ue) {
            response.setStatus(HttpServletResponse.SC_UNAUTHORIZED);
            response.addHeader("WWW-Authenticate", "Basic realm=\"Clotho Rest\"");
            response.addHeader("HTTP/1.0 401", "Unauthorized");
            response.getWriter().write("{\"error\": \"unauthorized access of page\"}");
            return;
        }

        String result = this.rc.getResult().toString();
        System.out.println("\n\n\n" + result + "\n\n\n");
        //logout(unamePass);

        response.getWriter().write(result);

        if (result.contains("FAILURE")) {
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
        } else {
            response.setStatus(HttpServletResponse.SC_OK);
        }
        
        logout(auth);
        
    }

    
    
    
    
    private void login(String[] userPass) {
        if (userPass != null) {
            loginMap = new HashMap<String, String>();
            loginMap.put("username", userPass[0]);
            loginMap.put("credentials", userPass[1]);
            loginMessage = new Message(Channel.login, loginMap, null, null);
            this.router.receiveMessage(this.rc, loginMessage);
        }
    }

    private void logout(String[] userPass) {
        if (userPass != null) {
            logoutMessage = new Message(Channel.logout, loginMap, null, null);
            this.router.receiveMessage(this.rc, logoutMessage);
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
        } catch (NullPointerException ne) {
            map = new HashMap<String, String>();
        }

        return map;
    }

    private Map<String, String> splitQuery(String queryInfo) {
        Map<String, String> queryPairs = new HashMap<String, String>();
        String[] pairs = queryInfo.split("&");
        for (String pair : pairs) {
            int i = pair.indexOf("=");
            queryPairs.put(pair.substring(0, i), pair.substring(i + 1));
        }
        return queryPairs;
    }
}

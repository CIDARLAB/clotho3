package org.clothocad.webserver.jetty;

import org.clothocad.core.communication.Channel;
import org.clothocad.core.communication.Message;
import org.clothocad.core.communication.RestConnection;
import org.clothocad.core.communication.Router;
import org.clothocad.core.persistence.Persistor;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import javax.inject.Inject;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import org.apache.shiro.authz.UnauthorizedException;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.util.JSON;
import org.clothocad.model.Person;
import org.clothocad.model.Sequence;
import org.json.JSONObject;
import org.apache.shiro.SecurityUtils;

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
    @Inject
    public RestApi(Persistor persistor, Router router) {
        this.router = router;
        this.persistor = persistor;
    }

    protected void doGet(HttpServletRequest request,
            HttpServletResponse response) throws ServletException, IOException {

        response.addHeader("Access-Control-Allow-Origin", "*");
        response.setContentType("application/json");

        String[] pathID = request.getPathInfo().split("/");
        String id = pathID[2];


        JSONObject body = getRequestBody(request.getReader());
        Collection<ObjBase> raw = persistor.listAll();

        String username = body.get("username").toString();
        String password = body.get("password").toString();
        String[] auth = {username, password};

        if(!login(auth)) 
        {
            response.getWriter().write("Login Failed\r\n");
            return;
        }

//        Map<String, Object> query = new HashMap<>();
//        query.put("name", "B0034 Sequence"); //List should include BBa_K249006
//        Iterable<ObjBase> rawtwo = persistor.findRegex(query);
//        for (ObjBase each : rawtwo) {
//            System.out.println("REGEX LIST : " + each);
//        }
        System.out.println("\n\n\n");

        //String data = pathID[3];
        switch (id) {
            case "autocomplete":
                break;

            case "startsWith":
                break;

            case "get":
//                String type = body.get("type");
                Map<String, Object> query = new HashMap<>();

//                switch (type) {
//                    case "sequence":
                String name = body.get("name").toString();
                query.put("name", name);

//                }
                Iterable<ObjBase> rawtwo = persistor.findRegex(query);

                for (ObjBase each : rawtwo) {
                    System.out.println("REGEX LIST : " + each);
                }

                break;

            case "getAll":
                break;

            case "learn":
                break;

            case "query":
                break;

            case "queryOne":
                break;

            case "validate":
                break;
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
        String id = pathID[2];


        JSONObject body = getRequestBody(request.getReader());
        Collection<ObjBase> raw = persistor.listAll();

        String username = body.get("username").toString();
        String password = body.get("password").toString();
        String[] auth = {username, password};

        if(!login(auth)) 
        {
            response.getWriter().write("Login Failed\r\n");
            return;
        }

        switch (id) {
            case "destroy":
                System.out.println("\n\n\n get \n\n\n");
                m = new Message(Channel.destroy, body, null, null);
                break;
            case "destroyAll":
                break;
            case "clear":
                break;
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

        response.getWriter().write(result);

        if (result.contains("FAILURE")) {
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
        } else {
            response.setStatus(HttpServletResponse.SC_OK);
        }

        logout(auth);
    }

    protected void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
        response.addHeader("Access-Control-Allow-Origin", "*");
        response.setContentType("application/json");

        String[] pathID = request.getPathInfo().split("/");
        String id = pathID[2];

        JSONObject body = getRequestBody(request.getReader());
        Collection<ObjBase> raw = persistor.listAll();

        String username = body.get("username").toString();
        String password = body.get("password").toString();
        String[] auth = {username, password};

        if (id.equals("createUser")) {
            Map<String, String> credentials = new HashMap<>();
            credentials.put("username", auth[0]);
            credentials.put("credentials", auth[1]);
            credentials.put("displayname", auth[0]);
            m = new Message(Channel.createUser, credentials, null, null);
            if(!login(auth)) 
            {
                response.getWriter().write("Login Failed\r\n");
                return;
            }
            logout(auth);
        }

        if(!login(auth)) 
        {
            response.getWriter().write("Login Failed\r\n");
            return;
        }

        // Elowitz RBS sequence
//        Sequence seqB0034 = new Sequence("B0034 Sequence", "aaagaggagaaa", user);
//        persistor.save(seqB0034);
//
//        System.out.println("\n\n\n Make \n\n\n");
//
//        for (ObjBase each : raw) {
//            System.out.println("ALL LIST : " + each);
//        }
//        System.out.println("\n\n\n");

        switch (id) {
            case "create":
                Person user = new Person(auth[0]);
                String type = body.get("type").toString();
                String name = body.get("name").toString();
                String value = body.get("value").toString();

                switch (type) {
                    case "sequence":
                        Sequence seq = new Sequence(name, value, user);
                        persistor.save(seq);
                        break;
                }

                break;
            case "createAll":
                m = new Message(Channel.createAll, body, null, null);
                break;
            case "convert":
                break;
            case "log":
                break;
            case "run":
                break;
            case "say":
                break;
            case "submit":
                break;
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
        String id = pathID[2];

        JSONObject body = getRequestBody(request.getReader());
        Collection<ObjBase> raw = persistor.listAll();

        String username = body.get("username").toString();
        String password = body.get("password").toString();
        String[] auth = {username, password};

        if(!login(auth)) 
        {
            response.getWriter().write("Login Failed\r\n");
            return;
        }

        switch (id) {
            case "changePassword":
                break;
            case "grant":
                break;
            case "set":
                break;
            case "setAll":
                break;
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

        response.getWriter().write(result);

        if (result.contains("FAILURE")) {
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
        } else {
            response.setStatus(HttpServletResponse.SC_OK);
        }

        logout(auth);
    }

    private boolean login(String[] userPass) {
        if (userPass != null) {
            loginMap = new HashMap<String, String>();
            loginMap.put("username", userPass[0]);
            loginMap.put("credentials", userPass[1]);
            loginMessage = new Message(Channel.login, loginMap, null, null);
            this.router.receiveMessage(this.rc, loginMessage);
        }
        return SecurityUtils.getSubject().isAuthenticated();
    }

    private void logout(String[] userPass) {
        if (userPass != null) {
            logoutMessage = new Message(Channel.logout, loginMap, null, null);
            this.router.receiveMessage(this.rc, logoutMessage);
        }
    }

    private JSONObject getRequestBody(BufferedReader reader) throws IOException {
        
        StringBuilder buffer = new StringBuilder();
        String line;
        while ((line = reader.readLine()) != null) {
            buffer.append(line);
        }
        String data = buffer.toString();
        
        return new JSONObject(data);
        
    }
}

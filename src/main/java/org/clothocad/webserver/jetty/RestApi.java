package org.clothocad.webserver.jetty;

import com.fasterxml.jackson.databind.ObjectMapper;
import org.clothocad.core.communication.Channel;
import org.clothocad.core.communication.Message;
import org.clothocad.core.communication.RestConnection;
import org.clothocad.core.communication.Router;
import org.clothocad.core.persistence.Persistor;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.HashSet;
import javax.inject.Inject;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import org.apache.shiro.authz.UnauthorizedException;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.util.JSON;
import org.clothocad.model.*;
import org.json.JSONObject;
import org.apache.shiro.SecurityUtils;
import org.clothocad.core.datums.ObjectId;
import static org.clothocad.webserver.jetty.ConvenienceMethods.createPart;
import org.json.JSONArray;

@SuppressWarnings("serial")
public class RestApi extends HttpServlet {

    /*
    This is history. Never delete this comment:
    
    REST API servlet is at idontremember.url/data/*    
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
        String method = pathID[2];

        JSONObject body = getRequestBody(request.getReader());
        Collection<ObjBase> raw = persistor.listAll();

        String username = body.get("username").toString();
        String password = body.get("password").toString();
        String[] auth = {username, password};

        if (!login(auth)) {
            response.getWriter().write("Login Failed\r\n");
            return;
        }

        String result = "";
        switch (method) {
            case "getByName":
                Map<String, Object> query = new HashMap<>();
                String objectName = body.get("objectName").toString();
                query.put("name", objectName);

                Iterable<ObjBase> rawtwo = persistor.find(query);
                ObjBase last = null;
                for (ObjBase each : rawtwo) {
                    last = each;
                }

                result = last.toString();

                break;

            case "getById":

                String id = body.getString("id");
                ObjectId objId = new ObjectId(id);
                Object obj = persistor.get(objId);
                result = obj.toString();
                break;
        }

        response.setStatus(HttpServletResponse.SC_OK);
        response.getWriter().write(result);

        if (result.contains("FAILURE")) {
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
        } else {
            response.setStatus(HttpServletResponse.SC_OK);
        }

        logout(auth);

    }

    protected void doDelete(HttpServletRequest request,
            HttpServletResponse response) throws ServletException, IOException {
        response.addHeader("Access-Control-Allow-Origin", "*");
        response.setContentType("application/json");

        String[] pathID = request.getPathInfo().split("/");
        String method = pathID[2];

        JSONObject body = getRequestBody(request.getReader());

        String username = body.getString("username");
        String password = body.getString("password");
        String[] auth = {username, password};

        if (!login(auth)) {
            response.getWriter().write("Login Failed\r\n");
            return;
        }

        String result = "";

        switch (method) {
            case "delete":
                if (persistor.has(new ObjectId(body.getString("id"))))
                {
                    persistor.delete(new ObjectId(body.getString("id")));
                    response.getWriter().write("Object has been deleted\r\n");
                }
                else response.getWriter().write("Object with id " + body.getString("id") + " does not exist\r\n");
                break;
        }

        response.setStatus(HttpServletResponse.SC_OK);
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
        String method = pathID[2];
        String type = pathID[3];

        JSONObject body = getRequestBody(request.getReader());
        Collection<ObjBase> raw = persistor.listAll();

        String username = body.getString("username");
        String password = body.getString("password");
        String[] auth = {username, password};

        if (method.equals("create") && type.equals("user")) {
            Map<String, String> credentials = new HashMap<>();
            credentials.put("username", auth[0]);
            credentials.put("credentials", auth[1]);
            credentials.put("displayname", auth[0]);
            m = new Message(Channel.createUser, credentials, null, null);
            this.router.receiveMessage(this.rc, m);
            if (!login(auth)) {
                response.getWriter().write("Login Failed\r\n");
                return;
            }
            logout(auth);
        }

        if (!login(auth)) {
            response.getWriter().write("Login Failed\r\n");
            return;
        }

        Person user = new Person(auth[0]);
        Sequence sequence = null;
        Part part = null;
        Feature feature = null;
        String objectName = "";
        String role = "";
        String sequenceName = "";

        String result = "";

        switch (method) {
            case "create":
                switch (type) {
                    case "sequence":
                        objectName = body.getString("objectName");
                        sequenceName = body.getString("sequence");
                        sequence = new Sequence(objectName, sequenceName, user);
                        ObjectId sequenceObj = persistor.save(sequence);
                        result = sequenceObj.toString();
                        break;

                    case "part":
                        objectName = body.getString("objectName");
                        if (body.has("id")) {
                            String sequenceId = body.getString("id");
                            ObjectId id = new ObjectId(sequenceId);
                            sequence = persistor.get(Sequence.class, id);
                            part = new Part(objectName, sequence, user);
                        } else {
                            part = new Part(objectName, user);
                        }

                        ObjectId partObj = persistor.save(part);
                        result = partObj.toString();
                        break;

                    case "feature":
                        objectName = body.getString("objectName");
                        role = body.getString("role");
                        feature = new Feature(objectName, role, user);
                        ObjectId featureObj = persistor.save(feature);
                        result = featureObj.toString();
                        break;

                    case "module":
                        objectName = body.getString("objectName");
                        role = body.getString("role");

                        if (body.has("id")) {
                            String featureId = body.getString("id");
                            ObjectId id = new ObjectId(featureId);
                            feature = persistor.get(Feature.class, id);
                        } else {
                            feature = new Feature(objectName, role, user);
                        }
                        Set<Feature> features = new HashSet<Feature>();
                        features.add(feature);
                        BasicModule module = new BasicModule(objectName, role, features, user);
                        ObjectId moduleObj = persistor.save(module);
                        result = moduleObj.toString();
                        break;
                }
                break;

            case "convenience":
                switch (type) {
                    case "createPart":
                        Map<String, String> params = new HashMap<>();
                        params.put("role", body.getString("role"));
                        params.put("sequence", body.getString("sequence"));
                        objectName = body.getString("objectName");

                        ObjectId partObj = createPart(persistor, objectName, params, user.toString());
                        result = partObj.toString();
                        break;
                }
                break;
        }

        response.setStatus(HttpServletResponse.SC_OK);
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
        String method = pathID[2];

        JSONObject body = getRequestBody(request.getReader());
        Collection<ObjBase> raw = persistor.listAll();

        String username = body.get("username").toString();
        String password = body.get("password").toString();
        String[] auth = {username, password};

        if (!login(auth)) {
            response.getWriter().write("Login Failed\r\n");
            return;
        }

        switch (method) {
            case "set":
                body.remove("username");
                body.remove("password");
                if (body.getString("id") == null) {
                    response.getWriter().write("You must supply an id \r\n");
                    break;
                }

                ObjectId id = new ObjectId(body.getString("id"));

                if (!persistor.has(id)) {
                    response.getWriter().write("No object with this id exists\r\n");
                    break;
                }

                persistor.save(jsonToMap(body));
                
                //Contact the user to notify them that they modified an object
                response.getWriter().write("Successfully modified object");

                break;
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
    
    public static Map<String, Object> jsonToMap(JSONObject json) {
        Map<String, Object> retMap = new HashMap<String, Object>();

        if(json != JSONObject.NULL) {
            retMap = toMap(json);
        }
        return retMap;
    }

    public static Map<String, Object> toMap(JSONObject object) {
        Map<String, Object> map = new HashMap<String, Object>();

        Iterator<String> keysItr = object.keys();
        while(keysItr.hasNext()) {
            String key = keysItr.next();
            Object value = object.get(key);

            if(value instanceof JSONArray) {
                value = toList((JSONArray) value);
            }

            else if(value instanceof JSONObject) {
                value = toMap((JSONObject) value);
            }
            map.put(key, value);
        }
        return map;
    }

    public static List<Object> toList(JSONArray array) {
        List<Object> list = new ArrayList<Object>();
        for(int i = 0; i < array.length(); i++) {
            Object value = array.get(i);
            if(value instanceof JSONArray) {
                value = toList((JSONArray) value);
            }

            else if(value instanceof JSONObject) {
                value = toMap((JSONObject) value);
            }
            list.add(value);
        }
        return list;
    }
}

package org.clothocad.webserver.jetty;

import org.clothocad.core.communication.Channel;
import org.clothocad.core.communication.Message;
import org.clothocad.core.communication.RestConnection;
import org.clothocad.core.communication.Router;
import org.clothocad.core.persistence.Persistor;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
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
import org.clothocad.core.datums.ObjBase;
import org.clothocad.model.*;
import org.json.JSONObject;
import org.apache.shiro.SecurityUtils;
import org.clothocad.core.datums.ObjectId;
import static org.clothocad.webserver.jetty.ConvenienceMethods.*;
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
        String toGet = pathID[3];
        String username = pathID[4].split(":")[0];
        String password = pathID[4].split(":")[1];

        String[] auth = {username, password};

//        if (!login(auth)) {
//            response.getWriter().write("Login Failed\r\n");
//            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
//            return;
//        }

        String result = "";
        switch (method) {
            case "getByName":
                Map<String, Object> query = new HashMap<>();
                query.put("name", toGet);

                Iterable<ObjBase> queried = persistor.find(query);
                ObjBase last = null;
                for (ObjBase each : queried) {
                    last = each;
                }

                if (last == null) {
                    response.setStatus(HttpServletResponse.SC_NOT_FOUND);
                } else {
                    result = last.toString();
                }

                break;

            case "getById":
                ObjectId objId = new ObjectId(toGet);
                Object obj = persistor.get(objId);

                if (obj == null) {
                    response.setStatus(HttpServletResponse.SC_NOT_FOUND);
                } else {
                    response.setStatus(HttpServletResponse.SC_FOUND);
                    result = obj.toString();
                }
                break;
        }

        response.getWriter().write(result);

//        logout(auth);
    }

    protected void doDelete(HttpServletRequest request,
            HttpServletResponse response) throws ServletException, IOException {
        response.addHeader("Access-Control-Allow-Origin", "*");
        response.setContentType("application/json");

        String[] pathID = request.getPathInfo().split("/");
        String method = pathID[2];
        String username = pathID[4].split(":")[0];
        String password = pathID[4].split(":")[1];

        JSONObject body = getRequestBody(request.getReader());

        String[] auth = {username, password};

//        if (!login(auth)) {
//            response.getWriter().write("Login Failed\r\n");
//            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
//            return;
//        }

        switch (method) {
            case "delete":
                if (persistor.has(new ObjectId(body.getString("id")))) {
                    persistor.delete(new ObjectId(body.getString("id")));
                    response.getWriter().write("Object has been deleted\r\n");
                    response.setStatus(HttpServletResponse.SC_OK);
                } else {
                    response.getWriter().write("Object with id " + body.getString("id") + " does not exist\r\n");
                    response.setStatus(HttpServletResponse.SC_NOT_FOUND);
                }
                break;
        }

//        logout(auth);
    }

    protected void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
        response.addHeader("Access-Control-Allow-Origin", "*");
        response.setContentType("application/json");

        String[] pathID = request.getPathInfo().split("/");
        String method = pathID[2];
        String type = pathID[3];
        String username = pathID[4].split(":")[0];
        String password = pathID[4].split(":")[1];

        JSONObject body = getRequestBody(request.getReader());

        String[] auth = {username, password};

        if (method.equals("create") && type.equals("user")) {
            Map<String, String> credentials = new HashMap<>();
            credentials.put("username", auth[0]);
            credentials.put("credentials", auth[1]);
            credentials.put("displayname", auth[0]);

            m = new Message(Channel.createUser, credentials, null, null);
            this.router.receiveMessage(this.rc, m);

//            if (!login(auth)) {
//                response.getWriter().write("Login Failed\r\n");
//                response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
//                return;
//            }
//            logout(auth);
        }

//        if (!login(auth)) {
//            response.getWriter().write("Login Failed\r\n");
//            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
//            return;
//        }

        Person user = new Person(auth[0]);
        Sequence sequence = null;
        Part part = null;
        Feature feature = null;
        String[] partIDs = {};
        String objectName = "";
        String role = "";
        String rawSequence = "";

        String result = "";

        switch (method) {
            case "create":
                switch (type) {
                    case "sequence":
                        objectName = body.getString("objectName");
                        rawSequence = body.getString("sequence");
                        sequence = new Sequence(objectName, rawSequence, user);
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

                    case "conveniencePart":
                        role = body.getString("role");
                        rawSequence = body.getString("sequence");
                        objectName = body.getString("objectName");

                        Map<String, String> partParams = new HashMap<>();
                        partParams.put("role", role);
                        partParams.put("sequence", rawSequence);

                        ObjectId partId = createPart(persistor, objectName, partParams, username);
                        result = partId.toString();
                        break;

                    case "convenienceDevice":
                        role = body.getString("role");
                        rawSequence = body.getString("sequence");
                        objectName = body.getString("objectName");
                        partIDs = body.getString("partIDs").split(",");
                        ArrayList<String> partIDArray = new ArrayList<>();

                        for (String partID : partIDs) {
                            partIDArray.add(partID);
                        }
                        
                        Map<String, String> deviceParams = new HashMap<>();
                        deviceParams.put("role", role);
                        deviceParams.put("sequence", rawSequence);

                        ObjectId deviceID = createDevice(persistor, objectName, partIDArray, deviceParams, username);
                        result = deviceID.toString();
                        break;
                }
                break;
        }

        if (result.isEmpty()) {
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
        } else {
            response.setStatus(HttpServletResponse.SC_CREATED);
        }

        response.getWriter().write(result);

//        logout(auth);
    }

    protected void doPut(HttpServletRequest request,
            HttpServletResponse response) throws ServletException, IOException {
        response.addHeader("Access-Control-Allow-Origin", "*");
        response.setContentType("application/json");

        String[] pathID = request.getPathInfo().split("/");
        String method = pathID[2];
        String username = pathID[4].split(":")[0];
        String password = pathID[4].split(":")[1];

        JSONObject body = getRequestBody(request.getReader());

        String[] auth = {username, password};

//        if (!login(auth)) {
//            response.getWriter().write("Login Failed\r\n");
//            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
//            return;
//        }

        switch (method) {
            case "set":
                body.remove("username");
                body.remove("password");
                if (!body.has("id")) {
                    response.getWriter().write("You must supply an id \r\n");
                    response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
                    break;
                }

                ObjectId id = new ObjectId(body.getString("id"));

                if (!persistor.has(id)) {
                    response.getWriter().write("No object with this id exists\r\n");
                    response.setStatus(HttpServletResponse.SC_NOT_FOUND);
                    break;
                }

                Map<String, Object> original = persistor.getAsJSON(id);

                Iterator<String> keysItr = body.keys();
                while (keysItr.hasNext()) {
                    String key = keysItr.next();
                    Object value = body.get(key);
                    original.put(key, value);
                }

                persistor.save(jsonToMap(body));
                response.setStatus(HttpServletResponse.SC_OK);

                break;
        }

//        logout(auth);
    }

//    private boolean login(String[] userPass) {
//        if (userPass != null) {
//            loginMap = new HashMap<String, String>();
//            loginMap.put("username", userPass[0]);
//            loginMap.put("credentials", userPass[1]);
//            loginMessage = new Message(Channel.login, loginMap, null, null);
//            this.router.receiveMessage(this.rc, loginMessage);
//        }
//        return SecurityUtils.getSubject().isAuthenticated();
//    }
//
//    private void logout(String[] userPass) {
//        if (userPass != null) {
//            logoutMessage = new Message(Channel.logout, loginMap, null, null);
//            this.router.receiveMessage(this.rc, logoutMessage);
//        }
//    }

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

        if (json != JSONObject.NULL) {
            retMap = toMap(json);
        }
        return retMap;
    }

    public static Map<String, Object> toMap(JSONObject object) {
        Map<String, Object> map = new HashMap<String, Object>();

        Iterator<String> keysItr = object.keys();
        while (keysItr.hasNext()) {
            String key = keysItr.next();
            Object value = object.get(key);

            if (value instanceof JSONArray) {
                value = toList((JSONArray) value);
            } else if (value instanceof JSONObject) {
                value = toMap((JSONObject) value);
            }
            map.put(key, value);
        }
        return map;
    }

    public static List<Object> toList(JSONArray array) {
        List<Object> list = new ArrayList<Object>();
        for (int i = 0; i < array.length(); i++) {
            Object value = array.get(i);
            if (value instanceof JSONArray) {
                value = toList((JSONArray) value);
            } else if (value instanceof JSONObject) {
                value = toMap((JSONObject) value);
            }
            list.add(value);
        }
        return list;
    }
}

package org.clothocad.webserver.jetty;

import com.google.common.collect.Lists;
import org.clothocad.core.communication.Channel;
import org.clothocad.core.communication.Message;
import org.clothocad.core.communication.RestConnection;
import org.clothocad.core.communication.Router;
import org.clothocad.core.persistence.Persistor;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
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
import org.clothocad.model.*;
import org.json.JSONObject;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.jongo.JongoConnection;
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
    private static Message m;
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

        String result = "";
        switch (method) {
            case "getByName":
                String lastId = pathID[4];
                String nextPrev = pathID[5];
                String pageSize = pathID[6];

                if (pageSize != null && !pageSize.isEmpty()) {

                    String query = "{name:\"" + toGet + "\"}";
                    String sortOrder = "";

                    //if going to the next page
                    if (nextPrev.equals("next") && (lastId != null && !lastId.isEmpty())) {
                        query = "{name:\"" + toGet + "\",_id:{ $gt : \"" + lastId + "\"}}";
                        sortOrder = "{_id:1}";
                    }

                    //if going to the previous page
                    if (nextPrev.equals("prev") && (lastId != null && !lastId.isEmpty())) {
                        query = "{name:\"" + toGet + "\",_id:{ $lte : \"" + lastId + "\"}}";
                        sortOrder = "{_id:-1}";
                    }
                    int per_page = Integer.parseInt(pageSize);
                    JongoConnection.Pagination queried = persistor.findByPage(query, sortOrder, per_page);

                    JSONObject jsono = new JSONObject();
                    jsono.put("page", queried.page);
                    jsono.put("per_page", queried.per_page);
                    jsono.put("page_count", queried.page_count);
                    jsono.put("total_count", queried.total_count);

                    JSONArray arr = new JSONArray();
                    JSONObject next = new JSONObject();

                    //if we are going to the previous page, the results will be in the wrong direction because of the way the query and sort order go
                    if (nextPrev.equals("prev") && (lastId != null && !lastId.isEmpty())) {
                        queried.list = Lists.reverse(queried.list);
                    }
                    //the last id in the records of the current page, used to navigate to next and previous pages
                    String newLastId = String.valueOf(queried.list.get(queried.list.size() - 1).getId());

                    next.put("next", "/" + newLastId + "/next/" + String.valueOf(pageSize));
                    if (queried.page != queried.page_count) {
                        arr.put(next);
                    }
                    JSONObject prev = new JSONObject();
                    prev.put("prev", "/" + newLastId + "/prev/" + String.valueOf(pageSize));
                    if (queried.page != 1) {
                        arr.put(prev);
                    }

                    jsono.put("links", arr);

                    jsono.put("records", new JSONArray(queried.list));

                    if (queried.list.isEmpty()) {
                        response.setStatus(HttpServletResponse.SC_NOT_FOUND);
                    } else {
                        response.setStatus(HttpServletResponse.SC_FOUND);
                        result = jsono.toString();
                    }
                } else {
                    System.out.println("Error: " + pageSize);
                    response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
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
    }

    protected void doDelete(HttpServletRequest request,
            HttpServletResponse response) throws ServletException, IOException {
        response.addHeader("Access-Control-Allow-Origin", "*");
        response.setContentType("application/json");

        String[] pathID = request.getPathInfo().split("/");
        String method = pathID[2];
        
        JSONObject body = getRequestBody(request.getReader());

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
    }

    protected void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
        response.addHeader("Access-Control-Allow-Origin", "*");
        response.setContentType("application/json");

        String[] pathID = request.getPathInfo().split("/");
        String method = pathID[2];
        String type = pathID[3];

        JSONObject body = getRequestBody(request.getReader());;

        if (method.equals("create") && type.equals("user")) {
            Map<String, String> credentials = new HashMap<>();
            credentials.put("username", body.get("username").toString());
            credentials.put("credentials", body.get("password").toString());
            credentials.put("displayname", body.get("username").toString());

            m = new Message(Channel.createUser, credentials, null, null);
            this.router.receiveMessage(this.rc, m);
        }

        String name = body.get("username").toString();
        Person user = new Person(name);
        Sequence sequence = null;
        Part part = null;
        Feature feature = null;
        String[] partIDs = {};
        String objectName = "";
        String role = "";
        String rawSequence = "";
        String result = "";
        JSONArray paramsArray = null;
        List params = null;
        Map<String, String> sequenceRole = null;


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

                        paramsArray = body.getJSONArray("params");
                        params = new ArrayList();
                        for (int i = 0; i < paramsArray.length(); i++) {
                            JSONObject childObject = paramsArray.getJSONObject(i);
                            String paramName = childObject.getString("name");
                            Double paramValue = childObject.getDouble("value");
                            String paramVariable = childObject.getString("variable");
                            String paramUnits = childObject.getString("units");
                            Parameter p = new Parameter(paramName, paramValue, paramVariable, paramUnits);
                            params.add(p);
                        }

                        sequenceRole = new HashMap<>();
                        sequenceRole.put("role", role);
                        sequenceRole.put("sequence", rawSequence);

                        ObjectId partId = createPart(persistor, objectName, sequenceRole, params, name);
                        result = partId.toString();
                        break;

                    case "convenienceDevice":
                        role = body.getString("role");
                        rawSequence = body.getString("sequence");
                        objectName = body.getString("objectName");
                        boolean createSeqFromParts = body.getBoolean("createSeqFromParts");
                        
                        paramsArray = body.getJSONArray("params");
                        params = new ArrayList();
                        for (int i = 0; i < paramsArray.length(); i++) {
                            JSONObject childObject = paramsArray.getJSONObject(i);
                            String paramName = childObject.getString("name");
                            Double paramValue = childObject.getDouble("value");
                            String paramVariable = childObject.getString("variable");
                            String paramUnits = childObject.getString("units");
                            Parameter p = new Parameter(paramName, paramValue, paramVariable, paramUnits);
                            params.add(p);
                        }
                        
                        
                        partIDs = body.getString("partIDs").split(",");
                        ArrayList<String> partIDArray = new ArrayList<>();
                        for (String partID : partIDs) {
                            partIDArray.add(partID);
                        }
                        
                        sequenceRole = new HashMap<>();
                        sequenceRole.put("role", role);
                        sequenceRole.put("sequence", rawSequence);

                        ObjectId deviceID = createDevice(persistor, objectName, partIDArray, sequenceRole, params, name, createSeqFromParts);
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
    }

    protected void doPut(HttpServletRequest request,
            HttpServletResponse response) throws ServletException, IOException {
        response.addHeader("Access-Control-Allow-Origin", "*");
        response.setContentType("application/json");

        String[] pathID = request.getPathInfo().split("/");
        String method = pathID[2];

        JSONObject body = getRequestBody(request.getReader());

        switch (method) {
            case "set":
                body.remove("username");
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

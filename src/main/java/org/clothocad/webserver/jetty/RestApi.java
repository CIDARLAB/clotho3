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
        int pathLength = pathID.length;

        String method = (pathLength >= 3) ? pathID[2] : null;
        String toGet = (pathLength >= 4) ? pathID[3] : null;
        String lastID = (pathLength >= 5) ? pathID[4] : null;
        String nextPrev = (pathLength >= 6) ? pathID[5] : null;
        String pageSize = (pathLength >= 7) ? pathID[6] : null;

        String result = "";
        switch (method) {
            case "getByName":
                if (pageSize != null && !pageSize.isEmpty()) {

                    String query = "{name:\"" + toGet + "\"}";
                    String sortOrder = "";

                    //if going to the next page
                    if (nextPrev.equals("next") && (lastID != null && !lastID.isEmpty())) {
                        query = "{name:\"" + toGet + "\",_id:{ $gt : \"" + lastID + "\"}}";
                        sortOrder = "{_id:1}";
                    }

                    //if going to the previous page
                    if (nextPrev.equals("prev") && (lastID != null && !lastID.isEmpty())) {
                        query = "{name:\"" + toGet + "\",_id:{ $lte : \"" + lastID + "\"}}";
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
                    if (nextPrev.equals("prev") && (lastID != null && !lastID.isEmpty())) {
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

            case "getByID":
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

        response.getWriter().write("\n" + result);
    }

    protected void doDelete(HttpServletRequest request,
            HttpServletResponse response) throws ServletException, IOException {
        response.addHeader("Access-Control-Allow-Origin", "*");
        response.setContentType("application/json");

        String[] pathID = request.getPathInfo().split("/");
        int pathLength = pathID.length;

        String method = (pathLength >= 3) ? pathID[2] : null;
        String id = (pathLength >= 4) ? pathID[3] : null;

        switch (method) {
            case "delete":
                if (persistor.has(new ObjectId(id))) {
                    persistor.delete(new ObjectId(id));
                    response.getWriter().write("\n Object has been deleted\r\n");
                    response.setStatus(HttpServletResponse.SC_OK);
                } else {
                    response.getWriter().write("\n Object with id " + id + " does not exist\r\n");
                    response.setStatus(HttpServletResponse.SC_NOT_FOUND);
                }
                break;
        }
    }

    protected void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
        response.addHeader("Access-Control-Allow-Origin", "*");
        response.setContentType("application/json");

        String[] pathID = request.getPathInfo().split("/");
        int pathLength = pathID.length;

        String method = (pathLength >= 3) ? pathID[2] : null;
        String type = (pathLength >= 4) ? pathID[3] : null;
        String option = (pathLength >= 5) ? pathID[4] : null;

        JSONObject body = getRequestBody(request.getReader());;

        // Input
        String username = body.has("username") ? body.getString("username") : null;
        String password = body.has("password") ? body.getString("password") : null;
        String objectName = body.has("objectName") ? body.getString("objectName") : null;
        String rawSequence = body.has("sequence") ? body.getString("sequence").toLowerCase() : null;
        String description = body.has("description") ? body.getString("description") : null;
        String role = body.has("role") ? body.getString("role") : null;
        String displayID = body.has("displayID") ? body.getString("displayID") : null;
        String[] partIDs = body.has("partIDs") ? body.getString("partIDs").split(",") : null;

        ObjectId id = body.has("id") ? new ObjectId(body.getString("id")) : null;

        JSONArray paramsArray = body.has("params") ? body.getJSONArray("params") : null;
        JSONArray partsArray = body.has("parts") ? body.getJSONArray("parts") : null;

        // Created
        String result = "";

        Sequence sequence = null;
        Person person = new Person(username);
        Part part = null;
        Feature feature = null;

        List params = null;
        List parts = null;

        Map<String, String> sequenceRole = null;
        Map<String, String> query = null;

        switch (method) {
            case "create":
                switch (type) {
//                  Needs error codes
                    case "user":
                        Map<String, String> credentials = new HashMap<>();
                        credentials.put("username", username);
                        credentials.put("credentials", password);
                        credentials.put("displayname", username);

                        m = new Message(Channel.createUser, credentials, null, null);
                        this.router.receiveMessage(this.rc, m);
                        result = "createUser";
                        break;

                    case "sequence":
                        sequence = new Sequence(objectName, description, rawSequence, person);
                        ObjectId sequenceObjID = persistor.save(sequence);
                        result = sequenceObjID.toString();
                        break;

                    case "part":
                        if (id != null) {
                            sequence = persistor.get(Sequence.class, id);
                            part = new Part(objectName, description, sequence, person);
                        } else {
                            part = new Part(objectName, description, person);
                        }

                        ObjectId partObj = persistor.save(part);
                        result = partObj.toString();
                        break;

                    case "feature":
                        feature = new Feature(objectName, description, role, person);
                        ObjectId featureObj = persistor.save(feature);
                        result = featureObj.toString();
                        break;

                    case "module":
                        if (id != null) {
                            feature = persistor.get(Feature.class, id);
                        } else {
                            feature = new Feature(objectName, role, person);
                        }
                        Set<Feature> features = new HashSet<Feature>();
                        features.add(feature);
                        BasicModule module = new BasicModule(objectName, description, role, features, person);
                        ObjectId moduleObj = persistor.save(module);
                        result = moduleObj.toString();
                        break;

                    case "conveniencePart":
                        if (paramsArray != null) {
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
                        }

                        if (role != null && rawSequence != null) {
                            sequenceRole = new HashMap<>();
                            sequenceRole.put("role", role);
                            sequenceRole.put("sequence", rawSequence);
                        }

                        if (paramsArray != null) {
                            if (role != null && rawSequence != null) {
                                result = createPart(persistor, objectName, displayID, sequenceRole, params, username).getValue();
                            } else {
                                result = createPart(persistor, objectName, displayID, params, username).getValue();
                            }
                        } else if (role != null && rawSequence != null) {
                            result = createPart(persistor, objectName, displayID, sequenceRole, username).getValue();
                        } else {
                            result = createPart(persistor, objectName, displayID, username).getValue();
                        }
                        break;

                    case "convenienceDevice":
                        boolean createSeqFromParts = body.getBoolean("createSeqFromParts");

                        ArrayList<String> partIDArray = new ArrayList<>();
                        for (String partID : partIDs) {
                            partIDArray.add(partID);
                        }

                        if (paramsArray != null) {
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
                        }

                        if (role != null && rawSequence != null) {
                            sequenceRole = new HashMap<>();
                            sequenceRole.put("role", role);
                            sequenceRole.put("sequence", rawSequence);
                        }

                        if (paramsArray != null) {
                            if (role != null && rawSequence != null) {
                                result = createDevice(persistor, objectName, displayID, partIDArray, sequenceRole, params, username, createSeqFromParts).getValue();
                            } else {
                                result = createDevice(persistor, objectName, displayID, partIDArray, params, username, createSeqFromParts).getValue();
                            }
                        } else if (role != null && rawSequence != null) {
                            result = createDevice(persistor, objectName, displayID, partIDArray, sequenceRole, username, createSeqFromParts).getValue();
                        } else {
                            result = createDevice(persistor, objectName, displayID, partIDArray, username, createSeqFromParts).getValue();
                        }
                        break;
                }
                break;

            case "get":

                Map<String, List> subObjects = new HashMap<String, List>();
                List<String> ids;
                JSONObject jsonDeviceRes;

                switch (type) {
                    case "conveniencePart":
                        Map<String, Map<String, String>> part_map = new HashMap<String, Map<String, String>>();

                        query = new HashMap<String, String>();
                        query.put("name", objectName);
                        query.put("displayID", displayID);
                        query.put("role", role);
                        query.put("sequence", rawSequence);

                        if (paramsArray != null) {
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
                        } else {
                            params = null;
                        }
                        option = (option != null && option.equalsIgnoreCase("exact")) ? "true" : "false";
                        part_map = getPart(persistor, query, params, Boolean.valueOf(option));

                        JSONObject jsonRes = new JSONObject(part_map);
                        result = jsonRes.toString();
                        break;

                    case "conveniencePartID":

                        query = new HashMap<String, String>();
                        query.put("name", objectName);
                        query.put("displayID", displayID);
                        query.put("role", role);
                        query.put("sequence", rawSequence);

                        if (paramsArray != null) {
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
                        } else {
                            params = null;
                        }
                        option = (option != null && option.equalsIgnoreCase("exact")) ? "true" : "false";
                        ids = getPartID(persistor, query, params, Boolean.valueOf(option));

                        jsonRes = new JSONObject(ids);
                        result = jsonRes.toString();
                        break;

                    case "convenienceDevice":

                        Map<String, Map<String, String>> device_map = new HashMap<String, Map<String, String>>();
                        query = new HashMap<String, String>();
                        query.put("name", objectName);
                        query.put("displayID", displayID);
                        query.put("role", role);
                        query.put("sequence", rawSequence);

                        if (paramsArray != null) {
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
                            subObjects.put("params", params);
                        }

                        if (partsArray != null) {
                            parts = new ArrayList();
                            for (int i = 0; i < partsArray.length(); i++) {
                                JSONObject childObject = partsArray.getJSONObject(i);
                                String partName = childObject.getString("name");
                                String partDescription = childObject.has("description") ? childObject.getString("description") : null;

                                Part p = new Part(partName, partDescription, person);
                                parts.add(p);
                            }
                            subObjects.put("parts", parts);
                        }

                        option = (option != null && option.equalsIgnoreCase("exact")) ? "true" : "false";
                        device_map = getDevice(persistor, query, subObjects, Boolean.valueOf(option));

                        jsonDeviceRes = new JSONObject(device_map);
                        result = jsonDeviceRes.toString();
                        break;

                    case "convenienceDeviceID":

                        query = new HashMap<String, String>();
                        query.put("name", objectName);
                        query.put("displayID", displayID);
                        query.put("role", role);
                        query.put("sequence", rawSequence);

                        if (paramsArray != null) {
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
                            subObjects.put("params", params);
                        }

                        if (partsArray != null) {
                            parts = new ArrayList();
                            for (int i = 0; i < partsArray.length(); i++) {
                                JSONObject childObject = partsArray.getJSONObject(i);
                                String partName = childObject.getString("name");
                                String partDescription = childObject.has("description") ? childObject.getString("description") : null;

                                Part p = new Part(partName, partDescription, person);
                                parts.add(p);
                            }
                            subObjects.put("parts", parts);
                        }
                        
                        option = (option != null && option.equalsIgnoreCase("exact")) ? "true" : "false";
                        ids = getDeviceID(persistor, query, subObjects, Boolean.valueOf(option));

                        jsonDeviceRes = new JSONObject(ids);
                        result = jsonDeviceRes.toString();
                        break;
                }
                break;
        }

        if (result.isEmpty()) {
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
        } else {
            response.setStatus(HttpServletResponse.SC_CREATED);
        }

        response.getWriter().write("\n" + result);
    }

    protected void doPut(HttpServletRequest request,
            HttpServletResponse response) throws ServletException, IOException {
        response.addHeader("Access-Control-Allow-Origin", "*");
        response.setContentType("application/json");

        String[] pathID = request.getPathInfo().split("/");
        int pathLength = pathID.length;
        String method = (pathLength >= 3) ? pathID[2] : null;

        JSONObject body = getRequestBody(request.getReader());
        ObjectId id = body.has("id") ? new ObjectId(body.getString("id")) : null;

        switch (method) {
            case "set":
                if (id == null) {
                    response.getWriter().write("You must supply an id \r\n");
                    response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
                    break;
                }

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

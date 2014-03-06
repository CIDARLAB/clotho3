package org.clothocad.core.communication;

import static org.clothocad.core.communication.Channel.autocompleteDetail;
import static org.clothocad.core.communication.Channel.create;
import static org.clothocad.core.communication.Channel.destroy;
import static org.clothocad.core.communication.Channel.log;
import static org.clothocad.core.communication.Channel.submit;

import java.io.IOException;
import java.lang.ref.WeakReference;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import javax.inject.Inject;
import javax.inject.Singleton;
import lombok.extern.slf4j.Slf4j;
import lombok.Getter;
import org.apache.shiro.SecurityUtils;
import org.apache.shiro.subject.Subject;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.Sharable;
import org.clothocad.core.execution.Mind;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.util.JSON;

@Slf4j
@Singleton
public class Router {

    @Getter
    protected Persistor persistor;
   
    @Inject
    public Router(Persistor persistor) {
        minds = new HashMap<>();
        this.persistor = persistor;
    }

    // send message    
    public void sendMessage(ClientConnection connection, Message message) {
        try {
            log.debug(JSON.serialize(message));
        } catch (IOException e) {
            log.debug("failed to serialize message: {}", message);
        }
        connection.send(message);
    }

    // receive message
    public void receiveMessage(ClientConnection connection, Message request) {

        //bind context to request
        Subject subject = SecurityUtils.getSubject();
        Mind mind;
        if (subject.isAuthenticated()){
            mind = getAuthenticatedMind(subject.getPrincipal().toString());
            mind.setConnection(connection);
        } else {
            mind = getMind(connection);
        }
        ServerSideAPI api = new ServerSideAPI(mind, persistor, this, request.getRequestId());


        Object data = request.getData();
        
        Object response = null;
        try {
            switch (request.getChannel()) {
                case autocomplete:
                    api.autocomplete(data.toString());
                    break;
                case autocompleteDetail:
                    api.autocompleteDetail(data.toString());
                    break;
                case submit:
                    response = api.submit(data.toString());
                    break;
                case clear:
                    api.clear();
                    break;
                case login:
                    
                    Map map = (Map) data;
                    
                    response = api.login(map.get("username").toString(), map.get("password").toString());
                    break;
                case logout:
                    String key = SecurityUtils.getSubject().getPrincipal().toString();                    
                    response = api.logout();
                    authenticatedMinds.remove(key);
                    break;
                case changePassword:
                    api.changePassword(data.toString());
                    break;
                case learn:
                    api.learn(data);
                    break;

                case log:
                    api.log(data.toString());
                    break;
                case say:
                    api.say(data);
                    break;
                case note:
                    api.note(data.toString());
                    break;
                case alert:
                    api.alert(data.toString());
                    break;

                case get:
                    response = api.get(data);
                    break;
                case set:
                    response = api.set(JSON.mappify(data));
                    if (response == null) response = Void.TYPE;
                    break;
                case create:
                    response = api.create(data);
                    break;
                case destroy:
                    response = api.destroy(data);
                    if (response == null) response = Void.TYPE;
                    break;
                case query:
                    response = api.query(JSON.mappify(data));
                    break;
                    
                //TODO:
                case getAll:
                    response = api.getAll((List) data);
                    break;
                case createAll:
                    response = api.createAll((List) data);
                    break;
                case destroyAll:
                    api.destroyAll((List) data);
                    response = Void.TYPE;
                    break;
                case setAll:
                    api.setAll((List) data);
                    response = Void.TYPE;
                    break;
                case queryOne:
                    response = api.queryOne(JSON.mappify(data));
                    break;

                case run:
                    response = api.run(data);
                    break;
                case listen:
                    api.listen(data.toString());
                    break;
                case unlisten:
                    api.unlisten(data.toString());
                    break;
                default:
                    log.warn("Unknown channel {}", request.getChannel());
                    break;
            }
            
            if (response == Void.TYPE){
                connection.deregister(request.getChannel(), request.getRequestId());
            }
            else {
                Message message = new Message(request.getChannel(), response, request.getRequestId());
                connection.send(message);
            }
            
        } catch (Exception e) {
            //TODO: message client with failure
            api.say(e.getMessage(), ServerSideAPI.Severity.FAILURE, request.getRequestId());
            log.error(e.getMessage(), e);
        }
    }

    //Start JCA's hack of a pubsub, to be replaced by Ernst
   /* void publish(Sharable object) {
     try {
     System.out.println("Ernst, this needs to be implemented.  Push object via pubsub.");
     String uuid = object.getUUID().toString();
     Map<String, Object> msg = persistor.
     HashSet<WeakReference<ClientConnection>> targets = pubsub.get(uuid);
     for(WeakReference<ClientConnection> wr : targets) {
     ClientConnection conn = wr.get();
     if(conn==null) {
     continue;
     }
     try {
     sendMessage(conn, msg);
     } catch(Exception err) { }
     }
                        
                
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        
     }*/
    void publish(Map object) {
        //doesn't do anything.
        //Is here because Stephanie thinks we should publish the JSON data, not the actual object
    }
    private HashMap<String, HashSet<WeakReference<ClientConnection>>> pubsub = new HashMap<>();

    void register(ClientConnection connection, Sharable object) {
        String uuid = object.getId().toString();
        HashSet<WeakReference<ClientConnection>> existing = pubsub.get(uuid);
        if (existing == null) {
            existing = new HashSet<>();
        }

        existing.add(new WeakReference<>(connection));
        pubsub.put(uuid, existing);
    }
    //End JCA's hack of a pubsub, to be replaced by Ernst
    
    private Mind getMind(ClientConnection connection) {
        String id = connection.getId();
        if (minds.containsKey(id)) {
            Mind mind = minds.get(id);
            if (mind.getConnection() != connection){
                //XXX: this is probably disasterous in some edge cases
                //because jetty preserves the sesson id across websocket close/open, need to check to see if the connection object in the mind is stale
                mind.setConnection(connection);
            }
            return mind;
        }

        Mind mind = new Mind();

        mind.setConnection(connection);
        
        minds.put(id, mind);
        return mind;
    }
    
    private Map<String, Mind> minds;
    private Map<String, Mind> authenticatedMinds = new HashMap<>();

    private Mind getAuthenticatedMind(String username)  {
        //XXX: this whole method is janky
        if (authenticatedMinds.containsKey(username)){
            return authenticatedMinds.get(username);
        }
        
        
        Map<String,Object> query = new HashMap();
        query.put("username", username);
        query.put("className", Mind.class.getCanonicalName());
        try{
            Iterable<ObjBase> minds = persistor.find(query);

        Mind mind;
        
        if (!minds.iterator().hasNext()){
            mind = new Mind();
            authenticatedMinds.put(username, mind);
        } else {
            mind = (Mind) minds.iterator().next();
        }
        
        return mind;
                }catch (Exception ex){
            ex.printStackTrace();
            throw ex;
        }
    }

}

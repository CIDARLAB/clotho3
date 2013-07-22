package org.clothocad.core.layers.communication;

import com.google.inject.Guice;
import com.google.inject.Injector;
import java.lang.ref.WeakReference;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import javax.inject.Inject;
import lombok.extern.slf4j.Slf4j;
import org.bson.types.ObjectId;
import org.clothocad.core.datums.Doo;
import org.clothocad.core.datums.Sharable;
import static org.clothocad.core.layers.communication.Channel.autocompleteDetail;
import static org.clothocad.core.layers.communication.Channel.create;
import static org.clothocad.core.layers.communication.Channel.destroy;
import static org.clothocad.core.layers.communication.Channel.log;
import static org.clothocad.core.layers.communication.Channel.submit;

import org.clothocad.core.layers.communication.connection.ClientConnection;
import org.clothocad.core.layers.communication.mind.Mind;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.persistence.mongodb.MongoDBModule;
import org.clothocad.core.util.JSON;

@Slf4j
public class Router {

    /**
     * DOUBLE-CHECKED LOCKING *
     */
    private static volatile Router router;

    public static Router get() {
        Router result = router;
        if (result == null) {
            synchronized (Router.class) {
                result = router;
                if (result == null) {
                    router = result = new Router();
                }
            }
        }
        return result;
    }
    private ExecutorService pool;
    
    @Inject
    private Persistor persistor;
   
    public Router() {

        // create the thread-pool of ssAPI objects
        this.pool = Executors.newFixedThreadPool(100);
        minds = new HashMap<>();
    }

    // send message    
    public void sendMessage(ClientConnection connection, Message message) {
        log.debug(JSON.serialize(message));
        connection.send(message);
    }

    // receive message
    public void receiveMessage(ClientConnection connection, Message request) {
        RouterDoo doo = new RouterDoo();

        //bind context to request
        Mind mind = getMind(connection);
        ServerSideAPI api = mind.getAPI();


        Channel chanEnum = request.channel;
        Object data = request.data;
        String id = request.requestId;
        try {
            switch (chanEnum) {
                case autocomplete:
                    api.autocomplete(data.toString());
                    break;
                case autocompleteDetail:
                    api.autocompleteDetail(data.toString());
                    break;
                case submit:
                    api.submit(data.toString());
                    break;
                case clear:
                    api.clear();
                    break;
                case login:
                    //TODO
                    api.login(null, null);
                    break;
                case logout:
                    api.logout();
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
                    api.get(data, id);
                    break;
                case set:
                    api.set(JSON.mappify(data), id);
                    break;
                case create:
                    api.create(data, id);
                    break;
                case destroy:
                    api.destroy(data, id);
                    break;
                case query:
                    api.query(JSON.mappify(data), id);
                    break;

                case run:
                    api.run(data, id);
                    break;
                case show:
                    api.show(data.toString(), null, null, null);
                    break;
                case startTrail:
                    api.startTrail(data.toString());
                    break;
                case edit:
                    api.edit(data.toString());
                    break;
                case listen:
                    api.listen(data.toString());
                    break;
                case unlisten:
                    api.unlisten(data.toString());
                    break;
                default:
                    break;
            }
        } catch (Exception e) {
            //TODO: message client with failure
            api.say(e.getMessage(), ServerSideAPI.Severity.FAILURE, id);
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
        String uuid = object.getUUID().toString();
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
            if (mind.getClientConnection() != connection){
                //XXX: this is probably disasterous in some edge cases
                //because jetty preserves the sesson id across websocket close/open, need to check to see if the connection object in the mind is stale
                mind.setClientConnection(connection);
            }
            return mind;
        }

        Mind mind = new Mind();

        ServerSideAPI api = new ServerSideAPI(mind, persistor);
        mind.setAPI(api);
        mind.setClientConnection(connection);
        
        minds.put(id, mind);
        return mind;
    }
    
    private Map<String, Mind> minds;

    /**
     * The Doo's that manage any Client-derived message. Doo's handle even key
     * commands to avoid synchronization issues.
     */
    public class RouterDoo extends Doo {

        public RouterDoo() {
            super(null, false);
        }
        Message message;
        ObjectId mindId;
    }
}

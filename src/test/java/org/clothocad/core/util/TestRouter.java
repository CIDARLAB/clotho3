/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import com.google.inject.Injector;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.inject.Inject;
import org.clothocad.core.communication.Channel;
import org.clothocad.core.communication.Message;
import org.clothocad.core.communication.Router;
import org.clothocad.core.communication.ClientConnection;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.Persistor;

/**
 *
 * @author spaige
 */
public class TestRouter extends Router {

    @Inject
    public TestRouter(Persistor persistor,Injector inject) {
        super(persistor,inject);
    }

    @Override
    public void receiveMessage(ClientConnection connection, Message request) {
        if (request.getChannel() == Channel.reloadModels){
            Map<String,Object> query = new HashMap<>();
            query.put("schema", "org.clothocad.core.schema.BuiltInSchema");
            List<Map<String,Object>> results = persistor.findAsBSON(query);
            List ids = new ArrayList();
            for (Map<String,Object> result : results){
                ids.add(new ObjectId(result.get("id").toString()));
            }
            
            results = persistor.findAsBSON(new HashMap());
            
            for (Map<String,Object> result : results){
                ObjectId id = new ObjectId(result.get("id").toString());
                if (ids.contains(id)){
                    continue;
                }
                persistor.delete(id);
            }
            
            TestUtils.setupTestData(persistor);
            connection.send(new Message(
                request.getChannel(),
                "",
                request.getRequestId(),
                null
            ));
        } else {
            super.receiveMessage(connection, request);
        }
    }
}

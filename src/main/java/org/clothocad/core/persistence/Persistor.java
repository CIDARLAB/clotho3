/*
Copyright (c) 2010 The Regents of the University of California.
All rights reserved.
Permission is hereby granted, without written agreement and without
license or royalty fees, to use, copy, modify, and distribute this
software and its documentation for any purpose, provided that the above
copyright notice and the following two paragraphs appear in all copies
of this software.

IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.

THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
ENHANCEMENTS, OR MODIFICATIONS.
 */

package org.clothocad.core.persistence;

import com.google.common.cache.Cache;
import com.mongodb.BasicDBObject;
import com.mongodb.DBObject;
import com.mongodb.util.JSON;
import java.net.UnknownHostException;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.inject.Inject;
import lombok.Delegate;
import org.bson.BSONObject;
import org.bson.types.ObjectId;
import org.clothocad.core.aspects.Aspect;
import org.clothocad.core.aspects.JSONSerializer;
import org.clothocad.core.persistence.ClothoConnection;
import org.clothocad.core.persistence.mongodb.MongoDBConnection;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;

/**
 * @author jcanderson
 * 
 * Manages writing/reading objects to/from the DB, and will eventually also coordinate that with caching
 * Right now is a very simple pass-through class.
 * 
 * queries do not check the cache
 * 
 * fooAsBSON does not check or populate the cache - maintaining a separate BSON cache would be complicating but might help performance-wise
 * 
 * should probably get the cache to cooperate with the object deserializer (or pass the cache in for query methods)
 * possible refactoring: separate deserializer from connection?
 * 
 * having multiple threads have access to the same object sounds like a recipe for disaster -
 * Persistor should hand out copies
 * 
 *   
 */

//TODO: thread safety
//TODO: check out date created/modified/accessed bugs
public class Persistor implements Aspect {
    
    @Delegate
    private ClothoConnection connection;
    
    @Delegate
    private JSONSerializer serializer;
    
    private Cache cache;
    
    public Persistor(ClothoConnection connection){
        this(connection, null);
    }
    
    @Inject
    public Persistor(ClothoConnection connection, JSONSerializer serializer){
        this.connection = connection;
        this.serializer = serializer;
        //cache = new  CacheBuilder<ObjectId, ObjBase>().
    }

    public static Persistor get() {
        return singleton;
    }
    
    /**
     * JCA:  Yes, I know this is ugly.  You can change to some dependency injection thing later.
     */
    private static  Persistor singleton;
    static {
        try {
            MongoDBConnection conn = new MongoDBConnection();
            conn.connect();
            singleton = new Persistor(conn);
        } catch (UnknownHostException ex) {
            ex.printStackTrace();
            System.exit(0);
        }
    }

    public void persistFeature(HashMap<String, Integer> StoreGrams, String feature) {
        System.out.println("Stephanie:  features should be stored in a separate spot than ObjBases.  They aren't UUID-based.  They are are feature-word based");
        try {
            //Convert the StoreGrams to a JSONArray in a JSONObject
            JSONObject bolus = new JSONObject();
            JSONArray data = new JSONArray();
            for(String key : StoreGrams.keySet()) {
                Integer count = StoreGrams.get(key);
                JSONObject item = new JSONObject();
                item.put("command", key);
                item.put("count", count);
                data.put(item);
            }
            bolus.put("data", data);
            bolus.put("feature", feature);

            DBObject bson = ( DBObject ) JSON.parse( bolus.toString() );

            //Query the database for an existing feature entry
            HashMap query = new HashMap();
            query.put("feature", feature);
            BSONObject result = connection.getOneAsBSON(query);
            
            ObjectId oid = null;
            if(result!=null) {
                oid = (ObjectId) result.get("_id");
            } else {
                oid = ObjectId.get();
            }
        
        
            //Install the OID object
            bson.put("_id", oid);
            
            //Save it to the database and return the uuid
            connection.save(new BasicDBObject(bson.toMap()));
            
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public HashMap<String, Integer> loadFeature(String feature) {
        try {
            System.out.println("Stephanie:  feature persistence needs to be separated from ObjBase persistence");

            //Query the database for this feature entry
            HashMap query = new HashMap();
            query.put("feature", feature);
            BSONObject result = connection.getOneAsBSON(query);


            HashMap<String, Integer> out = new HashMap<String, Integer>();
            if(result==null) {
                return out;
            }


            //Transfer the data to a well-typed Map

            //need to re-do the parsing out of the database entry
            String array = result.get("data").toString();
            JSONArray ary = new JSONArray(array);
            
            for(int i=0; i<ary.length(); i++) {
                JSONObject item = ary.getJSONObject(i);
                String command = item.getString("command");
                int count = item.getInt("count");
                out.put(command, count);
            }

            return out;
        } catch(Exception err) {
        }
        return null;
    }        
}

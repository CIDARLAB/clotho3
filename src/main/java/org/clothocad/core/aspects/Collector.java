/*
Copyright  (c) 2010 The Regents of the University of California.
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

package org.clothocad.core.aspects;

import org.clothocad.core.persistence.Persistor;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.bson.types.ObjectId;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.model.Person;
import org.json.JSONException;
import org.json.JSONObject;

/**
 * @author jcanderson
 */
public class Collector implements Aspect {

    /**
     * Getting objBases from cache
     */
    public ObjBase getObjBase (String uuid) {
        try {
            if (objBaseBag.containsKey(uuid)) {
                return objBaseBag.get(uuid);
            } else {
            	ObjBase obj = Persistor.get().get(ObjBase.class, new ObjectId(uuid));
            	return obj;
            }
        } catch(Exception err) {
            return null;
        }
    }

    /**
     * Log a ObjBase into the cache.
     * 
     * This also puts the ObjBase in the "oldest" list
     * @param obj 
     */
    public synchronized void add (ObjBase obj) {
        /* Logger.log(Logger.Level.INFO, obj.getUUID() + " is the obj"); */
        objBaseBag.put(obj.getUUID(), obj);
        dateAccessedMap.put(new Date(), obj);
        
        //Start dropping objbases if it's full
        if(objBaseBag.size()>MAXCACHESIZE || dateAccessedMap.size()>MAXCACHESIZE) {
            ObjBase oldestObjBase = dateAccessedMap.firstEntry().getValue();
            dateAccessedMap.remove(dateAccessedMap.firstEntry().getKey());
            objBaseBag.remove(oldestObjBase.getUUID());
        }
        /* TODO: shouldn't save each time, but rather flush periodically */
//        Persistor.get().persistObjBase(obj); //I don't know that this should be here
    }

    private Collector() { }
    
    public static Collector get() {
        return singleton;
    }

    private static final Collector singleton = new Collector();

    
    //this list of ObjBases currently in RAM
    private final HashMap<ObjectId,ObjBase> objBaseBag = new HashMap<ObjectId,ObjBase>();
    private final TreeMap<Date,ObjBase> dateAccessedMap = new TreeMap<Date,ObjBase>();

    //How many objBases to hold in RAM before dropping things
    private static final Integer MAXCACHESIZE = 20000;

    public ObjBase temporaryRefetchMethod(String uuid) {
        System.out.println("Stephanie should change the implemention here such that prexisting objects update their values by Reflection");
        ObjBase obj = Persistor.get().get(ObjBase.class, new ObjectId(uuid));
        objBaseBag.put(obj.getUUID(), obj);
        return obj;
    }
    
    public JSONObject getMetadata(String uuid) {
        try {
            JSONObject json = new JSONObject();
            json.put("executed",23);
            json.put("successful",223);
            json.put("positive",723);
            json.put("negative",3);
            return json;
        } catch (JSONException ex) {
            return null;
        }
    }
    public static Person getAdmin() {
            //If it is already puilled, return it
            if(admin!=null) {
                return admin;
            }
            
            //Try getting it from the db
            try {
                Map query = new HashMap();
                query.put("className", "org.clothocad.model.Person");
                query.put("name", "admin");
                List<ObjBase> listy = Persistor.get().get(query);
                admin = (Person) listy.get(0);
                if(admin!=null) {
                    return admin;
                }
            } catch (Exception ex) {
            }
        
            //Create a new person
            admin = new Person("admin", null, "");
            Persistor.get().save(admin);
            return admin;
    }
    
    
    private static transient Person admin; 
    
}

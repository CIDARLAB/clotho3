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

import java.util.Date;
import java.util.HashMap;
import java.util.TreeMap;

import org.clothocad.core.aspects.Aspect;
import org.clothocad.core.datums.Datum;
import org.clothocad.core.datums.objbases.Person;
import org.clothocad.core.util.Logger;

/**
 * @author jcanderson
 */
public class Collector implements Aspect {

    /**
     * Getting datums from cache
     */
    public Datum getDatum (String uuid) {
        try {
            if (datumBag.containsKey(uuid)) {
                return datumBag.get(uuid);
            } else {
                if(uuid.startsWith("static")) {                	
                    return getStaticDatum(uuid);
                } else {
                	return null;
                    //return Persistor.get().loadDatum(uuid);
                }
            }
        } catch(Exception err) {
            return null;
        }
    }

    private static Datum getStaticDatum(String uuid) {
        return staticDatums.get(uuid);
    }
    
    /**
     * Log a Datum into the cache.
     * 
     * This also puts the Datum in the "oldest" list
     * @param obj 
     */
    public synchronized void add (Datum obj) {
        /* Logger.log(Logger.Level.INFO, obj.getUUID() + " is the obj"); */
        datumBag.put(obj.getId(), obj);
        dateAccessedMap.put(new Date(), obj);
        
        //Start dropping objbases if it's full
        if(datumBag.size()>MAXCACHESIZE || dateAccessedMap.size()>MAXCACHESIZE) {
            Logger.log(Logger.Level.WARN, "Exceeded cache");
            Datum oldestDatum = dateAccessedMap.firstEntry().getValue();
            dateAccessedMap.remove(dateAccessedMap.firstEntry().getKey());
            datumBag.remove(oldestDatum.getId());
        }
        /* TODO: shouldn't save each time, but rather flush periodically */
//        Persistor.get().persistDatum(obj); //I don't know that this should be here
    }

    private Collector() { 
    }
    
    private static Collector collector = null;
    public static Collector get() {
    	if(null == collector) {
    		collector = new Collector();
    	}
    	return collector;
    }


    
    //this list of Datums currently in RAM
    private final HashMap<String,Datum> datumBag = new HashMap<String,Datum>();
    private final TreeMap<Date,Datum> dateAccessedMap = new TreeMap<Date,Datum>();

    //How many datums to hold in RAM before dropping things
    private static final Integer MAXCACHESIZE = 20000;
    
    //Pre-load the various static objects into a hashmap
    private static final HashMap<String, Datum> staticDatums = new HashMap<String, Datum>();
    static {
        staticDatums.put(Person.getSchema().getId(), Person.getSchema());
        staticDatums.put(Person.getAdmin().getId(), Person.getAdmin());
    }

}

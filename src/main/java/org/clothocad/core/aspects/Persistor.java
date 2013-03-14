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

package org.clothocad.core.aspects;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;

import org.bson.types.ObjectId;
import org.clothocad.core.aspects.Aspect;
import org.clothocad.core.datums.Datum;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.layers.persistence.flat.FlatFilePersistor;
import org.clothocad.core.layers.persistence.mongodb.MongoDBPersistor;

/**
 * @author jcanderson
 */
public abstract class Persistor 
	implements Aspect {

	private static final int N = 2;     
		// N==1 -> MongoDB persistor
		// N==2 -> FlotFile persistor
	
	public abstract boolean persist(Collection<ObjBase> col);	
    public abstract boolean persist(ObjBase obj);

    public abstract Datum get(Class c, ObjBase obj);
	public abstract Datum get(Class c, ObjectId id);
	public abstract ObjBase[] get(ObjBase obj);
	public abstract <T> T[] get(Class<T> t);
	public abstract boolean clearDB();	
	
	public abstract Datum loadDatum (String uuid);
    public abstract HashMap<String, Integer> loadFeature (String feature);
    public abstract List<String> loadWordBank ();
    public abstract void persistFeature (Object obj, String filePath);
    public abstract void persistWordBank (List<String> wordBank);
    
    
    public static Persistor get() {
    	switch(N) {
    	case 1:
    		return getMongoDBPersistor();    		
    	case 2:
    		return getFlatFilePersistor();
    	}
    	return (Persistor)null;
    }
    
    private static MongoDBPersistor mongo = null;
    private static MongoDBPersistor getMongoDBPersistor() {
    	if(null == mongo) {
    		mongo = new MongoDBPersistor();
    	}
    	return mongo;
    }

    private static FlatFilePersistor flat = null;
    private static FlatFilePersistor getFlatFilePersistor() {
    	if(null == flat) {
    		flat = new FlatFilePersistor();
    	}
    	return flat;
    }

}

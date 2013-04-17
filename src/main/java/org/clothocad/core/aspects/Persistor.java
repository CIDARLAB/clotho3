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

import java.net.UnknownHostException;
import java.util.Collection;
import javax.inject.Inject;
import lombok.Delegate;
import org.bson.types.ObjectId;
import org.clothocad.core.aspects.Aspect;
import org.clothocad.core.datums.Datum;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.layers.communication.Router;
import org.clothocad.core.layers.persistence.ClothoConnection;
import org.clothocad.core.layers.persistence.mongodb.MongoDBConnection;

/**
 * @author jcanderson
 * 
 * Manages writing/reading objects to/from the DB, and will eventually also coordinate that with caching
 * Right now is a very simple pass-through class.
 */
public class Persistor implements Aspect {
    
    @Inject
    public Persistor(ClothoConnection connection) throws UnknownHostException{
        this.connection = connection;
            connection.connect();
    }
    
    private interface Connect {
        void connect() throws UnknownHostException;
    }
    
    @Delegate(excludes=Connect.class)
    private ClothoConnection connection;
    
    
	private static volatile Persistor persistor;

    public static Persistor get() 
    		throws Exception {
    	Persistor result = persistor;
		if(result == null) {
			synchronized(Persistor.class) {
				result = persistor;
				if(result == null) {
					persistor = result = new Persistor(new MongoDBConnection());
				}
			}
		}
		return result;
    }
}

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
ENHANCEMENTS, OR MODIFICATIONS..
 */

package org.clothocad.core.layers.persistence;

import java.net.UnknownHostException;
import java.util.Collection;
import java.util.Date;
import java.util.List;
import org.bson.BSONObject;
import org.bson.types.ObjectId;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.util.ClothoDate;

/**
 *
 * @author J. Christopher Anderson
 * @author Douglas Densmore
 */
public interface ClothoConnection {

    /**
     * To tell the connection to connect to its database 
     * 
     * throws XXX if connection fails
     * 
     * @return
     */
    void connect() throws UnknownHostException;

    /**
     * Returns true if the database backing this connection adheres to the
     * Clotho data model.
     * XXX: what does it return if there is no connection?
     * @return
     */
    boolean isAClothoDatabase();

    /**
     * Programmatically tell the connection to disconnect itself
     * from its database
     * 
     * throws XXX if errors
     */
    void disconnect();

    /**
     * Determine whether the database is connected currently
     * @return true if currently connected
     */
    boolean isConnected();

    /**
     * Saves the given object and its references to the database.
     * 
     * throws XXX if errors
     * 
     * @param obj
     */
    void save(ObjBase obj);
    void save(BSONObject obj);
    
    /**
     * Saves the given collection of objects to the database.
     * @param objs
     * @return the number of objects successfully saved
     */
    int save(Collection<ObjBase> objs);

    /**
     * Delete the object from the database.
     * @param obj
     * @return
     */
    void delete( ObjBase obj );
    void delete(ObjectId id);

    /**
     * Deletes the given set of objects from the database.
     * @param objs
     * @return number of objects deleted
     */
    int delete( Collection<ObjBase> objs );
    

    /**
     * Returns the time the given ObjBase object was modified in the database.
     * @param obj
     * @return
     */
    ClothoDate getTimeModified( ObjBase obj );

    /**
     * Gets the object with the given uuid as the specified class,
     * or null if the object does not exist in the database.
     * 
     * @param type
     * @param uuid
     * @return
     */
    <T> T get(Class<T> type, ObjectId uuid);
    BSONObject getAsBSON(ObjectId uuid);

    /**
     * Search the database using a query using MongoDB semantics.
     * See: 
     * http://docs.mongodb.org/manual/core/read-operations/
     * http://docs.mongodb.org/manual/reference/operator/
     * 
     * Objects will be deserialized to their 'default' type (the type they were first saved as.)
     * 
     * Type information is stored in the 'type' field as the canonical java name of the default type
     * @param query
     * @return
     */
    Collection<ObjBase> get(BSONObject query );
    Collection<ObjBase> get(String name);
    Collection<BSONObject> getAsBSON(BSONObject query);
    Collection<BSONObject> getAsBSON(String name);    
    
    <T> T getOne(Class<T> type, BSONObject query);
    <T> T getOne(Class<T> type, String name);
    BSONObject getOneAsBSON(BSONObject query);
    BSONObject getOneAsBSON(String name);
    
    <T extends ObjBase> Collection<T> getAll(Class<T> type);
    
    //Deletes everything
    void deleteAll();
 
}

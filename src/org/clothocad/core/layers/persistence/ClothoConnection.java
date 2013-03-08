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

import java.util.Collection;
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
     * To tell the plugin to connect to its database, returns true
     * if it successfully connected
     * @return
     */
    boolean connect();

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
     * @return true if disconnected without errors
     */
    boolean disconnect();

    /**
     * Determine whether the database is connected currently
     * @return true if currently connected
     */
    boolean isConnected();

    /**
     * Saves the given object to the database.
     * Does not save object references - use ObjBase.save() for cascade save.
     * @param obj
     * @return true if the object was saved, false otherwise
     */
    boolean save( ObjBase obj );

    /**
     * Saves the given collection of objects to the database.
     * Does not save object references - use ObjBase.save() for cascade save.
     * @param objs
     * @return the number of objects successfully saved
     */
    int save( Collection<ObjBase> objs );

    // TODO: make delete callable only from the API
    /**
     * Delete the object from the database.
     * Should only be called by the API.
     * @param obj
     * @return
     */
    boolean delete( ObjBase obj );

    /**
     * Deletes the given set of objects from the database.
     * Should only be called by the API.
     * @param objs
     * @return
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
    <T> T get( Class<T> type, ObjectId uuid );
    BSONObject getAsBSON(ObjectId uuid);

    /**
     * Search the database using a query using MongoDB semantics.
     * 
     * Objects will be deserialized to their 'default' type (the type they were first saved as.)
     * 
     * Type information is stored in the 'type' field as the canonical java name of the default type
     * @param query
     * @return
     */
    Collection<ObjBase> get(BSONObject query );
    
    <T> T getOne(Class<T> type, BSONObject query);

    ObjBase[] getTableAsArray(ObjBase obj);
    
    boolean clear();

    /**
     * Returns an array of data for specificied types
     * Returns an array of arrays, one for each entry in the db of type type
     * The entries in the inner array are the uuid and name, and if the type is
     * Feature, its sequence, else an empty string
     * @param type
     * @return
     */
    <T> T[] getTableAsArray( Class<T> type );

    /**
     * Creates a criteria-based ClothoQuery object that can be used to
     * construct a specialized query to the database. Calling the getResults()
     * method on the ClothoQuery object will return the results of the query.
     * @param type
     * @return
     */
    ClothoQuery createQuery(Class type);

    //XXX: classless query
    
    /**
     * This interface represents the programmatic approach to creating complex
     * queries. It is based off Hibernate criteria based queries. The reason
     * for creating a new query system is to hide the hibernate layer from
     * the tools.
     *
     * All future connections which implement ClothoConnection should also
     * contain inner classes that implement ClothoQuery and extend
     * ClothoCriterion, respectively. ClothoQuery objects will be created
     * by the ClothoConnection, and outside code can then call the methods
     * described in this interface to configure the query before getting the
     * results using getResults.
     *
     * The query can be configured quickly using the methods that directly
     * return a ClothoQuery, such as between() and equal(). Every time
     * one of these methods is called, the new restriction should be added
     * as an add clause. For example, if we made the succesive calls
     * between( VOLUME, 5, 10 ) and equal( COLOR, "red" ), all resuls
     * returned should satisfy both of these properties.
     *
     * If you need to add a more complex criterion, you can create a
     * criterion using the ClothoCriterion class.
     */
    public interface ClothoQuery {

        /**
         * Return the type this query is performing its constraints on.
         * @return
         */
        public Class getType();

        /**
         * Creates a query that represents a constraint on the associated type
         * The current type of this query as returned by getType() must have a
         * reference to the argument type.
         * @param type
         * @return
         */
        public ClothoQuery createAssociationQuery( Enum field );

        /**
         * Sets the maximum number of results this query will return.
         * @param max
         * @return
         */
        public ClothoQuery setMaxResults( int max );

        /**
         * Adds the constraint that the values in the given field of the results
         * must be between hi and lo.
         * @param field
         * @param hi
         * @param lo
         * @return
         */
        public ClothoQuery between( Enum field, Object hi, Object lo );

        /**
         * Adds the constraint that the values in the given field of the results
         * must be less than hi.
         * @param field
         * @param hi
         * @return
         */
        public ClothoQuery lessThan( Enum field, Object hi );

        /**
         * Adds the constraint that the value in the given field of the results
         * must be greater than lo.
         * @param field
         * @param lo
         * @return
         */
        public ClothoQuery greaterThan( Enum field, Object lo );

        /**
         * Adds the constraint that the value in the given field of the results
         * must be equal to value.
         * @param field
         * @param value
         * @return
         */
        public ClothoQuery equal( Enum field, Object value );

        /**
         * Adds the constraint that the value in the given field of the results
         * must be equal to value.
         * @param field
         * @param value
         * @return
         */
        public ClothoQuery eq( Enum field, Object value );

        /**
         * Adds the constraint that the value in the given field of the results
         * must not be equal to value.
         * @param field
         * @param value
         * @return
         */
        public ClothoQuery notEqual( Enum field, Object value );

        /**
         * Adds the constraint that the value in the given field of the results
         * must not be equal to value.
         * @param field
         * @param value
         * @return
         */
        public ClothoQuery neq( Enum field, Object value );

        /**
         * Adds the constraint that the values in the given field of the results
         * must be null.
         * @param field
         * @return
         */
        public ClothoQuery isNull( Enum field );

        /**
         * Adds the constraint that the values in the given field of the results
         * must not be null.
         * @param field
         * @return
         */
        public ClothoQuery isNotNull( Enum field );

        /**
         * Adds the constraint that the values in the given field must match
         * the given pattern string. The pattern string can be a regex or SQL 
         * pattern string. For example, to search for values that contain
         * atcg, pass in "%atcg%" as value.
         * @param field
         * @param value
         * @return
         */
        public ClothoQuery matches( Enum field, String value, boolean caseSensitive );

        /**
         * Adds the constraint that the value in the given field must contain
         * the given string. 
         * @param field
         * @param value
         * @param caseSensitive
         * @return
         */
        public ClothoQuery contains( Enum field, String value, boolean caseSensitive );

        /**
         * Adds the given criterion to the set of restrictions on the
         * results returned by this query.
         * @param crit
         * @return
         */
        public ClothoQuery add( ClothoCriterion crit );

        /**
         * Adds the restriction that any of the given criterion must be
         * true for the results.
         * @param criteria
         * @return
         */
        public abstract ClothoCriterion or( Collection<ClothoCriterion> criteria );

        /**
         * Adds the restriction that either of the two criterion must
         * be true for the results.
         * @param a
         * @param b
         * @return
         */
        public abstract ClothoCriterion or(ClothoCriterion... criteria);

        /**
         * Adds the restriction that all of the givne criterion must be true
         * for the results.
         * @param criteria
         * @return
         */
        public abstract ClothoCriterion and( Collection<ClothoCriterion> criteria );

        /**
         * Adds the restriction that both of the given criterion must be
         * true for the results.
         * @param a
         * @param b
         * @return
         */
        public abstract ClothoCriterion and( ClothoCriterion... criteria );

        /**
         * Negates a ClothoCriterion
         * @param c
         * @return
         */
        public abstract ClothoCriterion not( ClothoCriterion c );

        /**
         * An object representing the condition that the value of the given
         * field must be between hi and lo. Can be used in the methods or(),
         * and(), and not() to create more complex ClothoCriterion, and then
         * added to a ClothoQuery using add().
         * @param field
         * @param hi
         * @param lo
         * @return
         */
        public abstract ClothoCriterion getBetweenCrit( Enum field, Object hi, Object lo );

        /**
         * An object representing the condition that the value of the given
         * field must be less than hi. Can be used in the methods or(),
         * and(), and not() to create more complex ClothoCriterion, and then
         * added to a ClothoQuery using add().
         * @param field
         * @param hi
         * @return
         */
        public abstract ClothoCriterion getLessThanCrit( Enum field, Object hi );

        /**
         * An object representing the condition that the value of the given
         * field must be greater than lo. Can be used in the methods or(),
         * and(), and not() to create more complex ClothoCriterion, and then
         * added to a ClothoQuery using add().
         * @param field
         * @param lo
         * @return
         */
        public abstract ClothoCriterion getGreaterThanCrit( Enum field, Object lo );

        /**
         * An object representing the condition that the value of the given
         * field must be equal to value. Can be used in the methods or(),
         * and(), and not() to create more complex ClothoCriterion, and then
         * added to a ClothoQuery using add().
         * @param field
         * @param value
         * @return
         */
        public abstract ClothoCriterion getEqualCrit( Enum field, Object value );

        /**
         * An object representing the condition that the value of the given
         * field must not be equal to value. Can be used in the methods or(),
         * and(), and not() to create more complex ClothoCriterion, and then
         * added to a ClothoQuery using add().
         * @param field
         * @param value
         * @return
         */
        public abstract ClothoCriterion getNotEqualCrit( Enum field, Object value );

        /**
         * An object representing the condition that the value of the given
         * field must be null. Can be used in the methods or(),
         * and(), and not() to create more complex ClothoCriterion, and then
         * added to a ClothoQuery using add().
         * @param field
         * @return
         */
        public abstract ClothoCriterion getIsNullCrit( Enum field );

        /**
         * An object representing the condition that the value of the given
         * field must not be null. Can be used in the methods or(),
         * and(), and not() to create more complex ClothoCriterion, and then
         * added to a ClothoQuery using add().
         * @param field
         * @return
         */
        public abstract ClothoCriterion getIsNotNullCrit( Enum field );

        /**
         * An object representing the condition that the value of the given
         * field must match value. Can be used in the methods or(),
         * and(), and not() to create more complex ClothoCriterion, and then
         * added to a ClothoQuery using add().
         * @param field
         * @param value
         * @return
         */
        public abstract ClothoCriterion getMatchesCrit( Enum field, Object value );

        /**
         * Once you are done creating the query, call this method to 
         * get the results. This method should not be called more than
         * once. 
         * @return
         */
        public List<Object> getResults();
    }

    public abstract class ClothoCriterion {
    }
}

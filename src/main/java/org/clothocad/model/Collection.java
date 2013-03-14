/*
Copyright (c) 2009 The Regents of the University of California.
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
package org.clothocad.model;

import com.github.jmkgreen.morphia.annotations.Reference;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.UUID;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.bson.types.ObjectId;
import org.clothocad.core.datums.ObjBase;

/**
 *
 * @author J. Christopher Anderson
 */
//name must be unique
public class Collection extends ObjBase {

    /**Constructor for collections from raw data
     *
     * @param collectionName = name of the Collection as a String
     * @param myauthor = author of Collection as a Person object
     */
    public Collection(String collectionName, String description, Person myauthor) {
        super(collectionName);
        this.description = description;
        this.author = myauthor;
        this.items = new ArrayList<ObjBase>();
    }

    /**Constructor for a transient Collection.  Transient collections
     * aren't saved to the database ever, but they can be passed around
     * via the Collector (they have UUIDs).
     *
     * @param collectionName = name of the Collection as a String
     * @param myauthor = author of Collection as a Person object
     */
    public Collection() {
        super("Results");
        description = "Transient Collection";
        items = new ArrayList<ObjBase>();
    }


    /**Abstract method addObject in general is for drag-and-drop events
     * Method gets called from the receiver of the drop
     *
     * @param dropObject is the object being dropped
     * @return is true if the drop is accepted for the receiver type
     */
    public boolean addObject(ObjBase dropObject) {
        if (items.contains(dropObject)) {
            return false;
        }

        AddAnyItem(dropObject);
        System.out.println("****Collection added an object " + dropObject.getName());
        return true;
    }

    /* SETTERS
     * */
    /**Remove an item from the Collection
     *
     * @param item = the item you want to remove
     */
    public boolean removeItem(ObjBase item) {
        if (!items.contains(item)){
            return false;
        }
        items.remove(item);
        return true;
    }

    /* GETTERS
     * */

    /**
     * Get all objects in the Collection
     * @return a List of ObjBase's in this Collection
     */
    public List<ObjBase> getAll() {
        List<ObjBase> out = new ArrayList<ObjBase>();
        for (ObjBase obj : items) {
            if (obj == null) {
                continue;
            }
            out.add(obj);
        }
        return out;
    }

    /**
     * Get all ObjBase in this Collection and any ObjBases in a Collection
     * within this Collection
     * @param type the type of ObjBase you want
     * @return a HashSet full of ObjBases
     */
    public <T extends ObjBase> Set<T> recursiveGetAllOf(Class<T> type) {
        return recursiveRelay(type, new HashSet<ObjectId>());
    }

    <T extends ObjBase> Set<T> recursiveRelay(Class<T> type, Set<ObjectId> tested) {
        tested.add(getUUID());
        Set<T> out = new HashSet<T>();
        List<Collection> allmycollections = getAll(Collection.class);
        for (Collection childcoll : allmycollections) {
            if (tested.contains(childcoll.getUUID())) {
                continue;
            }
            Set<T> itscontents = childcoll.recursiveRelay(type, tested);
            for (T o : itscontents) {
                out.add(o);
            }
        }

        for (T o : getAll(type)) {
            out.add(o);
        }

        return out;
    }


    /**
     * Get a list of UUIDs for ObjBase in this Collection of a given type
     * @param type the ObjType desired
     * @return a Set of ObjBases
     */
    public <T extends ObjBase> Set<UUID> getAllLinksOf(Class<T> type) {
        HashSet out = new HashSet<UUID>();
        for (ObjBase item : items){
            if (type.isInstance(item)) {
                out.add(item.getUUID());
            }
        }
        return out;
    }

    /**Generic getter for all of whatever type
     *
     * @param type the ObjType desired
     * @return a List of ObjBase of that type
     */
    public <T extends ObjBase> List<T> getAll(Class<T> type) {
        ArrayList out = new ArrayList<T>();
        for (ObjBase item : items){
            if (type.isInstance(item)) {
                out.add(item);
            }
        }
        return out;
    }

    /**
     * This one probably shouldn't exists
     * deprecated
     * @param _myPart
     * @return
     *
    Deprecated
    public ArrayList<Plasmid> getPlasmidsOf(Part _myPart) {
        @SuppressWarnings(value = "unchecked")
        ArrayList<Plasmid> allplas = (ArrayList<Plasmid>) getAll(ObjType.PLASMID);
        ArrayList<Plasmid> out = new ArrayList<Plasmid>();
        for (Plasmid p : allplas) {
            if (p.getPart().getUUID().equals(_myPart.getUUID())) {
                out.add(p);
            }
        }
        return out;
    }*/

    /**
     * This one probably shouldn't exists
     * deprecated
     * @param _myPart
     * @return
     *
    Deprecated
    public ArrayList<PlasmidSample> getSamplesOf(Plasmid _myPlasmid) {
        @SuppressWarnings(value = "unchecked")
        ArrayList<Sample> allsam = (ArrayList<Sample>) getAll(ObjType.SAMPLE);
        ArrayList<PlasmidSample> out = new ArrayList<PlasmidSample>();
        for (Sample p : allsam) {
            PlasmidSample ps = (PlasmidSample) p;
            System.out.println("comparing " + ps.getPlasmid().getUUID() + "  " + _myPlasmid.getUUID());

            if (ps.getPlasmid().getUUID().equals(_myPlasmid.getUUID())) {
                out.add(ps);
            }
        }
        return out;
    }*/

    private void AddAnyItem(ObjBase item) {
        if (!items.contains(item)){
            items.add(item);
        }
    }

    public void removeAll() {
        items.clear();
    }

    public static Collection retrieveByName(String name) {
        throw new UnsupportedOperationException();
    }


    /*-----------------
    variables
    -----------------*/
    
    @Getter
    @Setter        
    private String description;
    
    @Getter
    @Reference
    private Person author;
    @Reference
    private List<ObjBase> items;
            
   /* public static class CollectionDatum extends ObjBaseDatum {

        public Map<String, ObjType> uuidTypeHash = new HashMap<String, ObjType>();
        public Map<ObjType, HashSet<String>> typeUUIDHash = new EnumMap<ObjType, HashSet<String>>(ObjType.class);
        public ArrayList<String> itemUUIDs = new ArrayList<String>();
        public String _authorUUID;
        public String _description;

        @Override
        public ObjType getType() {
            return ObjType.COLLECTION;
        }
    }*/
}

package org.clothocad.core.layers.persistence.mongodb;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;

import org.bson.types.ObjectId;
import org.clothocad.core.aspects.Persistor;
import org.clothocad.core.datums.Datum;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.layers.persistence.ClothoConnection;

import com.mongodb.BasicDBObject;

public class MongoDBPersistor 
		extends Persistor {

    private ClothoConnection connection = null;

    public MongoDBPersistor() {
    	connection = new MongoDBConnection();
    	connection.connect();
    }
    
    /***
    //save children
    public boolean save() {
        return save(new HashSet<ObjBase>());
    }
    
    public boolean persist(Set<Datum> exclude){
    }
    
    //TODO: memoize (recursive version), 
    public static List<Field> getAllFields(Class c){
        ArrayList<Field> output = new ArrayList<Field>();
        while (c != null && c != Object.class){
            output.addAll(Arrays.asList(c.getDeclaredFields()));
            c = c.getSuperclass();
        }
        return output;
    }
    

    private static List<Field> getAllReferences(Class c){
        ArrayList<Field> output = new ArrayList<Field>();
        while (c != null && c != Object.class){
            for (Field f : c.getDeclaredFields()){
                if (f.getAnnotation(Reference.class) != null){
                    output.add(f);
                } 
            }
            c = c.getSuperclass();
        }
        return output;
    }
    
    public boolean delete() {
        return connection.delete(this);
    }
    
    //undo/redo commands
    //edit permissions
    
    //getters lazy eval references
    //setters update changed flag

	@Override
	public boolean persist(Collection<Datum> col) {
        boolean newId = this.UUID == null;
        if(newId) {
            this.UUID = ObjectId.get();
        }
        exclude.add(this);        
        ArrayList<ObjBase> children = new ArrayList<ObjBase>();
        
        for (Field f : getAllReferences(this.getClass())){
            boolean accessible = f.isAccessible();
            try {
                f.setAccessible(true);
                Object value = f.get(this);
                //reference might be a collection of references
                if (java.util.Collection.class.isInstance(value)){
                    //TODO: not typesafe
                    children.addAll((java.util.Collection) value);
                } else {
                    children.add((ObjBase) value);
                }
            } catch (IllegalArgumentException ex) {
                Logger.getLogger(ObjBase.class.getName()).log(Level.SEVERE, null, ex);
            } catch (IllegalAccessException ex) {
                Logger.getLogger(ObjBase.class.getName()).log(Level.SEVERE, null, ex);
            } finally {
                f.setAccessible(accessible);
            }
        }
        
        for (ObjBase child: children){
            //TODO: needs to detect and skip proxies if lazy loading is used
            if (!exclude.contains(child) && child != null){
                if (!child.save(exclude)){
                    return false;
                }
            }
        }
        
        //save self
        
        if  (!connection.save(this)){
            //if we set an id in anticipation of saving, but the save failed, revert to null (so id is consistent w/ db state)
            if (newId) this.UUID = null;
            return false;
        }
        return true;
	}
	**/
    
	@Override
	public boolean persist(Collection<ObjBase> col) {
		if(null != col) {
			for(ObjBase obj:col) {
				if(!this.persist(obj)) {
					return false;
				}
			}
			return true;
		}
		return false;
	}
    
	@Override
	public boolean persist(ObjBase obj) {
		this.connection.save(obj);
		return true;
	}

	@Override
	public Datum loadDatum(String uuid) {
		
		return null;
	}

	@Override
	public HashMap<String, Integer> loadFeature(String feature) {

		return null;
	}

	@Override
	public List<String> loadWordBank() {

		return null;
	}

	@Override
	public void persistFeature(Object obj, String filePath) {

	}

	@Override
	public void persistWordBank(List<String> wordBank) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public Datum get(Class c, ObjBase obj) {
        ObjectId id = obj.getUUID();
        if(null != id) {
            return this.get(c, id);
        } else {
            throw new UnsupportedOperationException("Not supported yet.");
        }        
	}

	public Datum get(Class c, ObjectId id) {
        return this.connection.get(c, id);
	}
	
	// this get() method returns an array of all T objects ...
	public ObjBase[] get(ObjBase obj) {		
		return this.connection.getTableAsArray(obj);
	}

	@Override
	public <T> T[] get(Class<T> t) {
		return this.connection.getTableAsArray(t);
	}

	@Override
	public boolean clearDB() {
		this.connection.clear();
		return false;
	}
}

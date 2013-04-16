
package org.clothocad.core.datums;

import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import lombok.AccessLevel;
import lombok.Data;
import lombok.NoArgsConstructor;
import lombok.Setter;

import org.bson.types.ObjectId;
import org.clothocad.core.aspects.Persistor;
import org.clothocad.core.datums.util.ClothoDate;
import org.json.JSONException;
import org.json.JSONObject;

import com.github.jmkgreen.morphia.annotations.Entity;
import com.github.jmkgreen.morphia.annotations.Id;
import com.github.jmkgreen.morphia.annotations.Reference;

/**
 *
 * @author spaige
 */
@Entity("data")
@Data
@NoArgsConstructor
public abstract class ObjBase 
		implements Datum {
	
	private static final long serialVersionUID = -1585051234016915083L;

	public ObjBase(String name) {
        this.name = name;
    }
    @Id
    private ObjectId UUID = null;
    
    //auto-populate dates for new entries without messing up deserialization
    //semantics - constructor for brand new entries, deserializer for preexisting entries
    //remove dateCreated getter
    @Setter(AccessLevel.PUBLIC)
    private String name;    
    private boolean isDeleted;    
	private ClothoDate dateCreated;
	private ClothoDate lastModified;
	private ClothoDate lastAccessed;
	
    protected JSONObject json;
	
	public ObjBase(JSONObject json) 
			throws JSONException {
		this.dateCreated = new ClothoDate();
		this.lastModified = new ClothoDate();
		this.lastAccessed = new ClothoDate();
		this.json = json;
		this.UUID = ObjectId.get();
	}
	
	/**
	public UUID getUUID() {
		return this.uuid;
	}
	**/
	
	public void update() {
		
		// here we need to call the client-side API
		// which forwards the update message 
		// to ``subscribed'' clients
		
	}
	
	public String getId() {
		if(null == this.UUID) {
			this.UUID = ObjectId.get();
		}
		return this.UUID.toString();
	}
	
	public ClothoDate getDateCreated() {
		return this.dateCreated;
	}
	
	public void setLastModified(ClothoDate lastModified) {
		this.lastModified = lastModified;
	}
	
	public ClothoDate getLastModified() {
		return this.lastModified;
	}
	
	
	public void setLastAccessed(ClothoDate lastAccessed) {
		this.lastAccessed = lastAccessed;
	}
	
	public ClothoDate getLastAccessed() {
		return this.lastAccessed;
	}
	
    // every Datum must implement the toJSON method
    public JSONObject toJSON()  {   
    	if (null != json) {
    		return this.json;
    	}
    	// here we can utilize Stephanie's introspection code to 
    	// map datums into JSON/BSON format
    	return new JSONObject();
    }
    
    public List<ObjBase> getChildren(){
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
        
        return children;
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
    
}


package org.clothocad.core.datums;

import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import lombok.AccessLevel;
import lombok.Data;
import lombok.NoArgsConstructor;
import lombok.Setter;

import org.bson.types.ObjectId;

import com.github.jmkgreen.morphia.annotations.Entity;
import com.github.jmkgreen.morphia.annotations.Id;
import com.github.jmkgreen.morphia.annotations.Reference;
import java.util.Date;
import java.util.Map;
import lombok.Getter;
import org.bson.BSONObject;
import org.clothocad.core.aspects.Persistor;
import org.json.JSONException;
import org.json.JSONObject;

/**
 *
 * @author spaige
 */
@Entity("data")
@Data
@NoArgsConstructor
public abstract class ObjBase {

    public ObjBase(String name) {
        this.name = name;
    }
    
    @Id
    private ObjectId UUID = new ObjectId();
    
    private String name;    
    private boolean isDeleted;    
    
    @Setter(AccessLevel.NONE)
    private Date dateCreated;
    private Date lastModified, lastAccessed;
	
	public void onUpdate() {
		// here we need to call the client-side API
		// which forwards the update message 
		// to ``subscribed'' clients
            
            //so, do all setters need to check to see if the value changed and then call onUpdate?
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
    
    public JSONObject toJSON() {
        
        //JCA's hack of re-pulling from db to serialize. Please change.
    	
    	// here, we should actually iterate over all data fields of the object
    	// (using java introspection for example)

    	System.out.println("[ObjBase.toJSON] -> "+this.UUID.toString());
        try {
            BSONObject bson = Persistor.get().getAsBSON(UUID);
            if(null != bson) {
            	return new JSONObject(bson.toString());
            }
            return new JSONObject();
        } catch (JSONException ex) {
            ex.printStackTrace();
            return null;
        }
    }
    
}

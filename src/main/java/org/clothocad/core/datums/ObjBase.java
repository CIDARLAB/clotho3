
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
import org.bson.BSONObject;
import org.clothocad.core.aspects.Collector;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.persistence.DBOnly;
import org.clothocad.model.Person;
import org.json.JSONObject;

/**
 *
 * @author spaige
 */
@Entity("data")
@Data
@NoArgsConstructor
public abstract class ObjBase implements Sharable {

    public ObjBase(String name) {
        this.name = name;
        
    }
    
    @Id
    private ObjectId UUID;
    private String name;
    private String author_id;
    
    @DBOnly
    private boolean isDeleted;    
    
    @Setter(AccessLevel.NONE)
    private Date dateCreated;
    @DBOnly
    private Date lastModified, lastAccessed;
	
	public void onUpdate() {
		
            
            
		// here we need to call the client-side API
		// which forwards the update message 
		// to ``subscribed'' clients
            
            //so, do all setters need to check to see if the value changed and then call onUpdate?
		
	}
        
    public List<ObjBase> getChildren(){
        ArrayList<ObjBase> children = new ArrayList<>();
        
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
    
    @Override
    public String getId() {
        return UUID.toString();
    }
    
    //TODO:
    @Override
    public JSONObject toJSON(){
        
        //JCA's hack of re-pulling from db to serialize.  Please change.
        try {
            //Pull the object from db, convert to JSONObject
            ObjectId uuid = this.getUUID();
            BSONObject bson = Persistor.get().getAsBSON(uuid);
            JSONObject out = new JSONObject(bson.toString());
            
            out.put("_id", "toberemoved");
            out.remove("_id");
            String uuidstr = uuid.toString();
            out.put("id", uuidstr);
            out.put("uuid", uuidstr); //this is a temporary hack to support Max's code.  It should just be id and a string.
            if(out.has("className")) {
                
/*  this is for how it should be:
   {
    "schema_id" : "asdgasdgt2q345",
   }
   * 
   * temporarily it has the className instead, but that's not correct
*/
                out.put("schema_id", out.getString("className"));
            }
            return out;
        } catch (Exception ex) {
            System.out.println("There appears to be some damaged data in your database, I'll ignore it");
            return null;
        }
    }

    @Override
    public Person getAuthor() {
        return (Person) Collector.get().getObjBase(author_id);
    }

    @Override
    /**
     * Override this in your ObjBase-extending class to customize the standard icon for this class of object
     */
    public String getIcon() {
        return "data:image/jpeg;base64,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";
    }

    @Override
    public String getDescfription() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
}

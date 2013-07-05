
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
import lombok.EqualsAndHashCode;
import org.clothocad.core.persistence.DBOnly;
import org.clothocad.core.persistence.Rename;

/**
 *
 * @author spaige
 */
@Entity("data")
@EqualsAndHashCode(exclude = {"dateCreated", "lastModified", "lastAccessed", "isDeleted"})
@Data()
@NoArgsConstructor
public abstract class ObjBase {
    
    //add schema
    //remove schema
    //can only manipulate schema set if you have write privs

    public ObjBase(String name) {
        this.name = name;
    }
    
    @Id
    @Rename("id")
    private ObjectId UUID;
    
    private String name;
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
        ArrayList<Field> output = new ArrayList<>();
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
    
    //TODO:
    /*private JSONObject toJSON(){
        
        //JCA's hack of re-pulling from db to serialize.  Please change.
        try {
            //Pull the object from db, convert to JSONObject
            ObjectId uuid = this.getUUID();
            BSONObject bson = persistor.getAsBSON(uuid);
            JSONObject out = new JSONObject(bson.toString());
            
            out.put("_id", "toberemoved");
            out.remove("_id");
            String uuidstr = uuid.toString();
            out.put("id", uuidstr);
            out.put("uuid", uuidstr); //this is a temporary hack to support Max's code.  It should just be id and a string.
            if(out.has("className")) {
    
/*  this is to deal with Max's funky:
    "$clotho" : {
        "schema" : "schema_person",
        "uuid" : "inst_second"
    },
*/
               // JSONObject dollarclotho = new JSONObject();  
               // dollarclotho.put("schema", out.getString("className"));
                //dollarclotho.put("uuid", uuidstr);
                //out.put("$clotho", dollarclotho);
                
/*  this is for how it should be:
   {
    "schema_id" : "org.clothod.models.Institution",
   }
*/ 
                /*
                out.put("schema_id", out.getString("className"));
            }
            return out;
        } catch (Exception ex) {
            System.out.println("There appears to be some damaged data in your database, I'll ignore it");
            return null;
        } 
    }*/
    
}

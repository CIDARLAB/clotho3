/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.schema;

import com.fasterxml.jackson.annotation.JsonView;
import java.util.Map;
import java.util.Set;
import javax.inject.Inject;
import lombok.Data;
import lombok.EqualsAndHashCode;
import lombok.NoArgsConstructor;
import lombok.extern.slf4j.Slf4j;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.persistence.annotations.Add;
import org.clothocad.core.persistence.annotations.Adds;
import org.clothocad.core.persistence.DBClassLoader;
import org.clothocad.core.persistence.jackson.JSONViews;
import org.clothocad.core.persistence.annotations.Reference;
import org.clothocad.model.Person;

/**
 *
 * @author spaige
 */

@EqualsAndHashCode(exclude={"fields", "methods"}, callSuper = true)
@Data
@NoArgsConstructor
@Slf4j
@Adds({@Add(name="language", provider="getLanguage"),
@Add(name="binaryName", provider="getBinaryName")})
public abstract class Schema extends SharableObjBase {
    
    public Schema(String name, String description, Person author){
        super(name, author, description);
    }
    
    protected static final String BASE_PACKAGE_BINARY = "org.clothocad.loadedschemas.";
   
    @Inject
    public static  DBClassLoader cl = null;
    
    @JsonView(JSONViews.Internal.class)
    protected byte[] classData;
    protected Map<String, ObjectId> dependencies;
    protected String source;
    
    //These are settable only in ClothoSchema - they are derived from source in other languages
    
    protected Set<ClothoField> fields;
    protected Set<Function> methods;
    
    @Reference
    protected Schema superClass;

    public abstract Language getLanguage();

    //TODO: handle files that result in multiple source files;
    //TODO: supply uuid's of referenced classes to resolve name conflicts
    public abstract void setSource(String source);
    //***Proxying methods to method handles???
    //can get bytecode from functions? 
   
    public String getBinaryName(){
        return BASE_PACKAGE_BINARY + "C"+ this.getId();
    }
    
    public static String getBinaryName(ObjectId id){
          return BASE_PACKAGE_BINARY + "C"+ id;      
    }
    
    public String getInternalName(){
        return getBinaryName().replace('.', '/');
    }
    
    //the classloader can only find saved schemas, so if this throws an exception, try saving the schema
    public Class<? extends ObjBase> getEnclosedClass(ClassLoader cl) throws ClassNotFoundException{
        return (Class<? extends ObjBase>) cl.loadClass(getBinaryName());
    }
    
    public static String extractIdFromClassName(String className){
        String[] a =  className.split("\\.");
        return a[a.length-1].substring(1);
    }
    
    public static boolean isSchemaClassName(String className){
        //XXX: needs to go away when schema class names become normal class names
        //Router router = Router.get();
        //if persistor is null, we're in bootstrapping, so don't bother
        /*if (router.getPersistor() == null) return false;
        if (router.getPersistor().resolveSchemaFromClassName(className) != null){
            return true;
        }*/
        return true;
    }  

    public boolean childOf(Schema schema) {
        if (schema == null || superClass == null) return false;
        if (schema.getId().equals(superClass.getId())) return true;
        return (childOf(schema.superClass));
    }
}

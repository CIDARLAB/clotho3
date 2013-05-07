/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.schema;

import java.util.Map;
import java.util.Set;
import lombok.Data;
import lombok.NoArgsConstructor;
import org.bson.types.ObjectId;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.Sharable;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.datums.util.Language;
import org.clothocad.model.Person;

/**
 *
 * @author spaige
 */
@Data
@NoArgsConstructor
public abstract class Schema extends Sharable {
    
    public Schema(String name, String description, Person author){
        super(name, author);
        this.description = description;
    }
    
    protected static final String BASE_PACKAGE_BINARY = "org.clothocad.loadedschemas.";
    
    protected byte[] classData;
    protected Map<String, ObjectId> dependencies;
    protected String description;
    protected String largeIconURL;
    protected String smallIconURL;
    protected String source;
    
    //These are settable only in ClothoSchema - they are derived from source in other languages
    
    protected Set<ClothoField> fields;
    protected Set<Function> methods;
    protected Schema superClass;

    
    public abstract Language getLanguage();

    //TODO: handle files that result in multiple source files;
    //TODO: supply uuid's of referenced classes to resolve name conflicts
    public abstract void setSource(String source);
    //***Proxying methods to method handles???
    //can get bytecode from functions? 
   
    public String getBinaryName(){
        return BASE_PACKAGE_BINARY + this.getUUID();
    }
    
    public String getInternalName(){
        return getBinaryName().replace('.', '/');
    }
    
    //the classloader can only find saved schemas, so if this throws an exception, try saving the schema
    public <T extends ObjBase> Class<T> getEnclosedClass(ClassLoader cl) throws ClassNotFoundException{
        return (Class<T>) cl.loadClass(getBinaryName());
    }
    
    public static String extractIdFromClassName(String className){
        String[] a =  className.split("\\.");
        return a[a.length-1];
    }
}

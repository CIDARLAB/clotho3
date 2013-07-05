/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.schema;

import com.github.jmkgreen.morphia.annotations.Reference;
import java.util.Map;
import java.util.Set;
import javax.inject.Inject;
import lombok.Data;
import lombok.EqualsAndHashCode;
import lombok.NoArgsConstructor;
import lombok.extern.slf4j.Slf4j;
import org.bson.types.ObjectId;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.Sharable;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.persistence.Add;
import org.clothocad.core.persistence.Adds;
import org.clothocad.core.persistence.DBClassLoader;
import org.clothocad.core.persistence.DBOnly;
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
public abstract class Schema extends Sharable {
    
    public Schema(String name, String description, Person author){
        super(name, author);
        this.description = description;
    }
    
    protected static final String BASE_PACKAGE_BINARY = "org.clothocad.loadedschemas.";
   
    @Inject
    public static  DBClassLoader cl = null;
    
    @DBOnly
    protected byte[] classData;
    protected Map<String, ObjectId> dependencies;
    protected String description;
    protected String largeIconURL;
    protected String smallIconURL;
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
        return BASE_PACKAGE_BINARY + "C"+ this.getUUID();
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
        return ObjectId.isValid(extractIdFromClassName(className));
    }  
}

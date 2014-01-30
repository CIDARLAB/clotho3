/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.datums;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import lombok.Getter;
import lombok.Setter;
import org.bson.types.ObjectId;
import org.clothocad.core.datums.util.Language;
import static org.clothocad.core.datums.util.Language.JAVASCRIPT;
import org.clothocad.core.execution.JavaScriptScript;
import org.clothocad.core.execution.Script;
import org.clothocad.core.persistence.Reference;

/**
 *
 * @author spaige
 */
public class Module extends ObjBase {
    
    protected Script code;
    @Reference
    protected Module[] dependencies;
    protected String description;
    @Getter
    @Setter
    protected Language language;

    public Module() {
    }
    
    
    public String encodeScript(){
        return code.toString();
    }
    
    public void setCode(String code){
        setCode(code, language);
    }
    
    protected void setCode(String code, Language language){
        this.code = generateScript(code,language);
    }
    
    public static Script generateScript(String code, Language language){
        switch (language){
            case JAVASCRIPT:
                return new JavaScriptScript(code);
            default:
                throw new UnsupportedOperationException("unsupported language");
        }          
    }
    
    public void decodeScript(Map obj){
        setCode((String) obj.get("code"), Language.valueOf((String) obj.get("language")));
    }
    
//    @PrePersist
//    protected void prePersist() {
//        syncDependencies();
//        //store code as plain string instead of script
//    }
//    
//    @PreLoad
//    protected void postLoad() {
//        
//    }
    
    //TODO: test w/ @PostLoad also
    protected void syncDependencies(){
        //figure out dependencies declared in code 
        Set<ObjectId> declaredDependencies = code.findDependencies();
        Set<ObjectId> listedDependencies = getDependencySet(dependencies);
      
        Set<ObjectId> listedButNotDeclared = new HashSet<>(listedDependencies);
        listedButNotDeclared.removeAll(declaredDependencies);
        //Set<ObjectId> declaredButNotListed;
        
        //declare all listed dependencies in code
        //code.addImports(listedButNotDeclared);
        //TODO: add declared but not listed dependencies to dependency list
    }
    
    private static Set<ObjectId> getDependencySet(Module[] dependencies) {
        Set<ObjectId> output = new HashSet<>();
        for (Module obj : dependencies){
            output.add(obj.getId());
        }
        return output;
    }
    
    public String getCode() {
        return code.getSource();
    }
    
    public String getCodeToLoad() {
        return code.encapsulateModule(code.getSource(), getSetup());
    }
    
    public String getSetup(){
        return this.code.generateImports(getDependencySet(dependencies));
    }

    public Function getFunction(String name) {
        Function function = new Function(name, null, null, null, language);
        function.dependencies = new Module[]{this};
        String formatString = "function () { return %s.%s.apply(this,arguments);};";
        function.setCode(String.format(formatString, this.getName(), name));
        return function;
    }
}

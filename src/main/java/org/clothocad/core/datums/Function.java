/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.datums;

import com.github.jmkgreen.morphia.annotations.Reference;
import javax.script.ScriptException;
import lombok.Getter;
import org.clothocad.core.datums.util.Language;
import org.clothocad.model.Person;

/**
 *
 * @author spaige
 */

//variable number of type parameterization?
//java object for first-class functions?
//declaration of requrirements?
public class Function extends ObjBase {
    
    public Function(){};
    
    public Function(String name, String[] argNames, Class[] input, Class[] output, String source, Language language){
        this.inputTypes = input;
        this.outputTypes = output;
        this.setName(name);
        //this.action = new ScriptEngineScript(name, source, language);
        this.argNames = argNames;
    }
    
    @Reference
    private Person author;
    private String description;
    
    @Getter
    private String[] argNames;
    
    private String[] dependencies;
    
    private Class[] inputTypes;
    //XXX: losing some duck-typing style flexibility
    private Class[] outputTypes;
    //XXX: if single return type, could make things more typesafe java-side
    //XXX: I don't even know what to do with this
    //XXX: all our target languages have single return value, so multiple return value is undefined
        //XXX: python has automatic tuple unpacking, is that what is intended?
    
    private Script preconditions;
    private Script action;
    
    public boolean canDooIt(Object... args){
        //args match input types
        try {
            for (int i = 0; i<inputTypes.length; i++){
                if (!inputTypes[i].isInstance(args[i])) return false;
                //XXX: throw type exception
             }
        } catch (IndexOutOfBoundsException e){
            return false;
        }

        
        if (preconditions != null) try {
            return (Boolean) preconditions.run(args);
        } catch (ScriptException ex) {
            return false;
        }
        return true;
    }
    
    public Object execute(Object... args) throws ScriptException{
        if (canDooIt(args)){
            return action.run(args);
        }
        return null;
    }
}

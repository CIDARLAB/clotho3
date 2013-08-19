/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.datums;

import com.github.jmkgreen.morphia.annotations.Reference;
import java.util.List;
import javax.script.ScriptException;
import lombok.Getter;
import lombok.Setter;
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
    
    public Function(String name, Argument[] arguments, Class output, String source, Language language){
        this.outputType = output;
        this.setName(name);
        //this.action = new ScriptEngineScript(name, source, language);
        this.args = arguments;
    }
    
    @Getter
    @Setter
    private Language language;
    
    @Reference
    private Person author;
    private String description;
    
    @Getter 
    private Argument[] args;
    
    @Getter 
    private FunctionTest[] tests;
    
    @Reference
    private Function[] dependencies;
    
    //XXX: losing some duck-typing style flexibility
    private Class outputType;
    //XXX: if single return type, could make things more typesafe java-side
    //XXX: I don't even know what to do with this
    //XXX: all our target languages have single return value, so multiple return value is undefined
        //XXX: python has automatic tuple unpacking, is that what is intended?
    
    private Script precondition;
    private Script code;
    
    public boolean canDooIt(Object... args){
        //args match input types
        try {
            for (int i = 0; i<this.args.length; i++){
                if (!this.args[i].getType().isInstance(args[i])) return false;
                //XXX: throw type exception
             }
        } catch (IndexOutOfBoundsException e){
            return false;
        }

        
        if (precondition != null) try {
            return (Boolean) precondition.run(args);
        } catch (ScriptException ex) {
            return false;
        }
        return true;
    }
    
    //TODO: convert to dict-style
    public Object execute(Object... args) throws ScriptException{
        if (canDooIt(args)){
            return code.run(args);
        }
        return null;
    }

    public static class FunctionTest {
        private  List<Object> args;
        private  Object value;

        public FunctionTest(List<Object> argValues, Object expectedResult) {
            this.args = argValues; 
            this.value = expectedResult;
        }
        
        public FunctionTest(){};
    }
    
    
}

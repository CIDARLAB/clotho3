/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.layers.execution;

import com.github.jmkgreen.morphia.annotations.Transient;
import com.google.common.base.Joiner;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.script.Invocable;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;

/**
 *
 * @author spaige
 */
public class JavaScriptScript implements Script{
    //needs arg names, function name
    
    //how well does one engine for all functions scale?
    //do we need to manage the number of global functions? 
    
    //TODO: each script function needs to execute in its own scope
    //but
    //script functions should be cached somehow
    public static ScriptEngine engine = ClothoScriptEngineManager.getEngineByName("javascript");
    
    public JavaScriptScript(){};
    
    public JavaScriptScript(String name, String source){
        this.source = source;
        this.name = name;
    }
    
    private String source;
    private String name;
    
    @Transient
    private boolean loaded = false;
    
    private void load() throws ScriptException{
        Object result = null;
        try {
            result = engine.get(name);
        } catch(IllegalArgumentException e){
        }
        if (result == null){
            engine.eval(source);
        }
        
        loaded = true;
    }
    
    public Object run(Object... args) throws ScriptException {
        Invocable invocable = (Invocable) engine;
        if (!loaded) load();
        try {
            Object function = engine.get(name);
            return invocable.invokeFunction(name, args);
        } catch (NoSuchMethodException ex) {
            //should never happen, thanks to load?
            return null;
        }
    }
    

    
}

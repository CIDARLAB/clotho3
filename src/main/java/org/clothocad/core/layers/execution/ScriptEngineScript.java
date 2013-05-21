/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.layers.execution;

import com.github.jmkgreen.morphia.annotations.Transient;
import com.github.jmkgreen.morphia.annotations.PostLoad;
import javax.script.Invocable;
import javax.script.ScriptEngine;
import javax.script.ScriptException;
import org.clothocad.core.datums.util.Language;

/**
 *
 * @author spaige
 */

public class ScriptEngineScript implements Script {
        //needs arg names, function name
    
    //how well does one engine for all functions scale?
    //do we need to manage the number of global functions? 
    
    //TODO: each script function needs to execute in its own scope
    //but
    //script functions should be cached somehow
    @Transient
    protected ScriptEngine engine;
       
    protected String source;
    protected String name;
    private Language language;
    
    public ScriptEngineScript(){};
    
    public ScriptEngineScript(String name, String source, Language language){
        this.name = name;
        this.source = source;
        this.language = language;
        
        getScriptEngine();
    }
    
    @Transient
    protected boolean loaded = false;
    
    protected void load() throws ScriptException{
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
    
    @PostLoad
    public void getScriptEngine(){
        engine = ClothoScriptEngineManager.getEngineByName(language.toString());
    }
    
    
    
    public Object run(Object... args) throws ScriptException {
        Invocable invocable = (Invocable) engine;
        if (!loaded) load();
        try {
            return invocable.invokeFunction(name, args);
        } catch (NoSuchMethodException ex) {
            throw new RuntimeException(ex);
        }
    }
}

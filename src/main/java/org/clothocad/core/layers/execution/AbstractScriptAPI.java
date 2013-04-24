/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.layers.execution;

import javax.script.ScriptEngine;
import org.bson.BSONObject;
import org.bson.BasicBSONObject;
import org.clothocad.core.datums.Function;
import org.clothocad.core.layers.persistence.ClothoConnection;


abstract public class AbstractScriptAPI implements ScriptAPI {
    
    
    @Override
    public void setEngine(ScriptEngine engine){
        this.engine = engine;
    }
    
    protected ScriptEngine engine;
    protected ClothoConnection connection;
    

    protected Function getFunctionByName(String functionName){
        BSONObject query = new BasicBSONObject("name", functionName);
        query.put("className","org.clothocad.core.execution.Function");
        
        return connection.getOne(Function.class, query);
    }
    
    //should check if already loaded a la Function
    public void importFunction(String name){
        Object result = null;
        try {
            result = engine.get(name);
        } catch(IllegalArgumentException e){
        }
        if (result == null){
            injectFunction(getFunctionByName(name));
        }
    }

    abstract public void injectFunction(Function f);
}

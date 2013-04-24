/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.layers.execution;

import javax.script.ScriptException;
import org.clothocad.core.datums.Function;


public class PythonAPI extends AbstractScriptAPI {

    public void injectFunction(Function f) {
        try {
            engine.put(f.getName(), f);
            engine.eval(f.getName() + " = " + f.getName() + ".execute");
        } catch (ScriptException ex) {
            throw new RuntimeException(ex);
        }
        
    }

    
}

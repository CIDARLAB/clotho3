/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.layers.execution;

import javax.script.ScriptEngine;

/**
 *
 * @author spaige
 */
public interface ScriptAPI {
    public void importFunction(String name);
    public void setEngine(ScriptEngine engine);
}

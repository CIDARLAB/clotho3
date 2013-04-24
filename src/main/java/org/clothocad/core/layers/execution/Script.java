/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.layers.execution;

import javax.script.ScriptException;

/**
 *
 * @author spaige
 */
public interface Script{
    public Object run(Object... args) throws ScriptException;

}

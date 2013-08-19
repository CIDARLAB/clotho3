/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.schema;

import java.lang.reflect.Method;

/**
 *
 * @author spaige
 */
public class ReflectionUtils {
    public static Method findMethodNamed(String name, int argsLength, Class type){
        for (Method method : type.getMethods()){
            if (method.getName().equals(name) && method.getParameterTypes().length == argsLength){
                return method;
            }
        }
        return null;
    }
}

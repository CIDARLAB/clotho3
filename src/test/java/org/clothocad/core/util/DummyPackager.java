/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import java.util.List;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.Module;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.datums.util.Language;

/**
 *
 * @author spaige
 */
public class DummyPackager extends Function {
    
    public DummyPackager(){
    }
    
    public DummyPackager(String name){
        setName(name);
        setId(new ObjectId(name));
        dependencies = new Module[0];
        description = "A module packager";
        language = Language.JAVA;    
    }
    
    public static DummyPackager createDummyPackager(){
        return new DummyPackager("packager");
    }
    
    public Module packJavaScript(String name, String description, String code, List<Module> dependencies){
        Module result;
        result = new Module(name, description, Language.JAVASCRIPT, code,
                dependencies.toArray(new Module[0]));   
        return result;
    }
}

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.execution;

import java.util.Collection;
import java.util.Set;
import lombok.Getter;
import org.clothocad.core.datums.ObjectId;

/**
 *
 * @author spaige
 */
public class PythonScript implements Script{
    
    public PythonScript(String name, String source){
        
    }

    public Object run(Object... args) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Getter
    private String source;

    @Override
    public Set<ObjectId> findDependencies() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public String generateImports(Collection<ObjectId> listedButNotDeclared) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public String modularizeFunction(String code) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public String encapsulateModule(String code, String setupCode) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }


    
}

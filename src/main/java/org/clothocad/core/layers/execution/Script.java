/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.layers.execution;

import java.util.Collection;
import java.util.Set;
import org.bson.types.ObjectId;

/**
 *
 * @author spaige
 */
public interface Script{
    public Set<ObjectId> findDependencies();
    
    public String getSource();

    public String generateImports(Collection<ObjectId> listedButNotDeclared);

}

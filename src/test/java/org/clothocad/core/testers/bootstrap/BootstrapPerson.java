/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.testers.bootstrap;

import org.clothocad.core.persistence.Persistor;
import org.clothocad.model.*;

/**
 * 
 * @author jcanderson
 */
public class BootstrapPerson {
    public static void main(String[] args) {
        Institution inst = new Institution("Brigham Young", "Utah City", "Utah", "USA");
        Persistor.get().save(inst);
        
        String id = inst.getId();
        System.out.println(id);
    }
}

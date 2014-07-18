/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import org.clothocad.core.datums.Module;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.datums.util.Language;

/**
 *
 * @author spaige
 */
public class SecurityTester extends Module {
    
    public SecurityTester(){
        setName(SecurityTester.class.getSimpleName());
        dependencies = new Module[0];
        language = Language.JAVA;
    }
    
    public String function(){
        return "function ran!";
    }
}

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.testers;

import java.util.Properties;
import org.clothocad.core.ClothoModule;

/**
 *
 * @author spaige
 */
public class ClothoTestModule extends ClothoModule {
    public ClothoTestModule(Properties props){
        super(props);
        defaults.put("dbname", "testClotho");
    }
    
    public ClothoTestModule(){
        this(null);
    }
}

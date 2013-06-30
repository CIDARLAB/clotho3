/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.testers;

import java.util.Properties;
import org.clothocad.core.persistence.mongodb.MongoDBModule;

/**
 *
 * @author spaige
 */
public class MongoDBTestModule extends MongoDBModule{
     public MongoDBTestModule(Properties props){
        super(props);
        defaults.put("dbName", "testClotho");
    }
    
    public MongoDBTestModule(){
        this(null);
    }   
}

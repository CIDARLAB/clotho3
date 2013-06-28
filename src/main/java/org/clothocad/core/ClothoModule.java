/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core;

import com.google.inject.AbstractModule;
import com.google.inject.name.Names;
import java.util.Properties;

/**
 *
 * @author spaige
 */
public class ClothoModule extends AbstractModule {
    public ClothoModule(Properties properties){
        this.properties = new Properties(defaults);
        this.properties.putAll(properties);
    }
    
    private static final Properties defaults = new Properties();
    static {
        defaults.put("port", "8080");
        defaults.put("dbname", "clotho");
    }
    
    private Properties properties;
    

    @Override
    protected void configure() {
        Names.bindProperties(binder(), properties);
    }
    
    
}

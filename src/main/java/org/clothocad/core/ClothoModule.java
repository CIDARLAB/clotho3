/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core;

import com.google.inject.AbstractModule;
import com.google.inject.name.Names;
import java.util.Properties;
import org.clothocad.core.layers.communication.Router;

/**
 *
 * @author spaige
 */
public class ClothoModule extends AbstractModule {
    public ClothoModule(Properties properties){
        defaults.put("port", "8080");
        defaults.put("dbname", "clotho");
        defaults.put("dbhost", "localhost");
        defaults.put("dbport", "27017");
        defaults.put("loglevel", "warn"); //TODO: doesn't affect anything yet
        this.properties = new Properties(defaults);
        if (properties != null) this.properties.putAll(properties);
    }
    
    public ClothoModule(){
        this(null);
    }
    
    protected final Properties defaults = new Properties();
    
    private final Properties properties;
    

    @Override
    protected void configure() {
        //TODO: make router completely DI (?)
        requestInjection(Router.get());
        Names.bindProperties(binder(), properties);
    }
    
    
}

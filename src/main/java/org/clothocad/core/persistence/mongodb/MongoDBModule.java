/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.mongodb;

import com.github.jmkgreen.morphia.mapping.Mapper;
import com.google.inject.AbstractModule;
import com.google.inject.Singleton;
import com.google.inject.name.Names;
import java.util.Properties;
import org.clothocad.core.aspects.JSONSerializer;
import org.clothocad.core.persistence.ClothoConnection;

/**
 *
 * @author spaige
 */
public class MongoDBModule extends AbstractModule {
    
    @Override
    protected void configure() {
        
        //set up morphia slf4j logging
        
        bind(ClothoConnection.class).to(MongoDBConnection.class);
        bind(JSONSerializer.class).to(ClothoMapper.class);
        bind(Mapper.class).to(ClothoMapper.class);
        bind(ClothoMapper.class).in(Singleton.class);
    }
}

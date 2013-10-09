/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core;

import com.google.inject.AbstractModule;
import com.google.inject.Provides;
import com.google.inject.name.Names;
import java.io.File;
import java.io.FileInputStream;
import java.security.KeyStore;
import java.util.Properties;
import org.clothocad.core.persistence.IdUtils;
import org.clothocad.core.schema.Schema;
import org.clothocad.core.security.SecurityModule;
import org.eclipse.jetty.servlet.ServletContextHandler;

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
        defaults.put("keystorepath", System.getProperty("java.home") + "/lib/security/cacerts".replace('/', File.separatorChar)); //;
        defaults.put("keystorepass", ""); //XXX: storing passwords as strings in memory is insecure
        defaults.put("keymanagerpass", "");
        this.properties = new Properties(defaults);
        if (properties != null) this.properties.putAll(properties);
    }
    
    public ClothoModule(){
        this(null);
    }
    
    protected final Properties defaults = new Properties();
    
    protected final Properties properties;
    

    @Override
    protected void configure() {
        Names.bindProperties(binder(), properties);
        requestStaticInjection(Schema.class);
        requestStaticInjection(IdUtils.class);
        
        //XXX: put this somewhere more reasonable
        ServletContextHandler sch = new ServletContextHandler();
        bind(ServletContextHandler.class).annotatedWith(Names.named("containerServletContext")).toInstance(sch);
        binder().install(new SecurityModule(sch));
    }
    
    @Provides
    protected KeyStore provideKeyStore() throws Exception {
        //loads keystore from provided parameters
        KeyStore keystore = KeyStore.getInstance(KeyStore.getDefaultType());
        char[] password = properties.getProperty("keystorepass").toCharArray();
        keystore.load(new FileInputStream(properties.getProperty("keystorepath")), password);
        return keystore;
    }
    
}

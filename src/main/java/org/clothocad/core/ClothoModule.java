/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core;

import com.google.inject.AbstractModule;
import com.google.inject.Provides;
import com.google.inject.name.Names;
import java.io.File;
import java.util.Properties;
import org.clothocad.core.persistence.IdUtils;
import org.clothocad.core.schema.Schema;
import org.clothocad.core.security.SecurityModule;
import org.eclipse.jetty.server.ssl.SslConnector;
import org.eclipse.jetty.server.ssl.SslSelectChannelConnector;
import org.eclipse.jetty.servlet.ServletContextHandler;
import org.eclipse.jetty.util.ssl.SslContextFactory;

/**
 *
 * @author spaige
 */
public class ClothoModule extends AbstractModule {
    public ClothoModule(Properties properties){
        for (ConfigOption option : ConfigOption.values()){
            defaults.put(option.name(), option.defaultValue);
        }
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
    protected SslConnector provideSslConnector() throws Exception {
        SslContextFactory cf = new SslContextFactory();
        cf.setKeyStorePath(properties.getProperty("keystorepath"));
        cf.setKeyStorePassword(properties.getProperty("keystorepass"));
        SslSelectChannelConnector sslConnector = new SslSelectChannelConnector(cf);                
        return sslConnector;
    }
    
}

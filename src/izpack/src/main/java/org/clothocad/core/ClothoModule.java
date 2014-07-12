/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core;

import com.google.inject.AbstractModule;
import com.google.inject.Provides;
import com.google.inject.name.Names;
import java.util.Properties;
import javax.inject.Singleton;
import org.clothocad.core.persistence.DBClassLoader;
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
    private final Properties config;

    public ClothoModule(Properties config) {
        this.config = config == null ?
            ConfigOption.getDefaultConfig() : config;
    }

    @Override
    protected void configure() {
        Names.bindProperties(binder(), config);
        requestStaticInjection(Schema.class);
        requestStaticInjection(IdUtils.class);

        bind(DBClassLoader.class).in(Singleton.class);
        
        //XXX: put this somewhere more reasonable
        ServletContextHandler sch = new ServletContextHandler();
        bind(ServletContextHandler.class).annotatedWith(Names.named("containerServletContext")).toInstance(sch);
        binder().install(new SecurityModule(sch));
    }
    
    @Provides 
    protected SslConnector provideSslConnector() throws Exception {
        SslContextFactory cf = new SslContextFactory();
        cf.setKeyStorePath(config.getProperty("keystorepath"));
        cf.setKeyStorePassword(config.getProperty("keystorepass"));
        SslSelectChannelConnector sslConnector = new SslSelectChannelConnector(cf);                
        return sslConnector;
    }
    
}

package org.clothocad.core;

import org.clothocad.core.persistence.DBClassLoader;
import org.clothocad.core.persistence.IdUtils;
import org.clothocad.core.schema.Schema;

import com.google.inject.AbstractModule;
import com.google.inject.name.Names;
import com.google.inject.Provides;

import org.eclipse.jetty.server.ssl.SslConnector;
import org.eclipse.jetty.server.ssl.SslSelectChannelConnector;
import org.eclipse.jetty.util.ssl.SslContextFactory;

import java.util.Properties;

import javax.inject.Singleton;

/**
 * General-purpose Guice module.
 * Processes config options and handles singletons
 */
public class ClothoModule extends AbstractModule {

    protected final Properties config;

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
    }

    /**
     * Provides an SslConnector backed by a KeyStore to Guice.
     * Uses the 'keystorepath' and 'keystorepass' config options
     */
    @Provides
    protected SslConnector provideSslConnector() throws Exception {
        SslContextFactory cf = new SslContextFactory();
        cf.setKeyStorePath(config.getProperty("keystorepath"));
        cf.setKeyStorePassword(config.getProperty("keystorepass"));
        SslSelectChannelConnector sslConnector = new SslSelectChannelConnector(cf);
        return sslConnector;
    }
}

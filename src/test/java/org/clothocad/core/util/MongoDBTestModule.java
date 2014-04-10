/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import com.google.inject.Singleton;
import org.clothocad.core.persistence.ClothoConnection;
import org.clothocad.core.persistence.mongodb.MongoDBModule;
import org.clothocad.core.security.CredentialStore;

/**
 *
 * @author spaige
 */
public class MongoDBTestModule extends MongoDBModule {

    @Override
    protected void configure() {
        super.configure(); //To change body of generated methods, choose Tools | Templates.
        bind(TestEnvConnection.class).in(Singleton.class);
    }

    
    
    @Override
    protected void configureConnection() {
        bind(ClothoConnection.class).to(TestEnvConnection.class);
    }

    @Override
    protected void configureCredentialStore() {
        bind(CredentialStore.class).to(TestEnvConnection.class);
    }
    
}

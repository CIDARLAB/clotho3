/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.utils;

import org.clothocad.core.security.SecurityModule;
import com.google.inject.Guice;
import com.google.inject.Injector;
import com.google.inject.servlet.ServletModule;
import java.util.Properties;
import org.apache.shiro.SecurityUtils;
import org.apache.shiro.mgt.SecurityManager;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.persistence.mongodb.MongoDBModule;
import org.clothocad.core.testers.ClothoTestModule;
import org.clothocad.webserver.jetty.ClothoWebserver;

/**
 * Starts a clotho server instance configured for testing
 * @author spaige
 */
public class ClothoTestEnvironment {
        public static void main(String[] args)
            throws Exception {

        Integer nPort = 8080;
        if (args.length > 1) {
            nPort = Integer.parseInt(args[0]);
        }
        
        Properties properties = new Properties();
        properties.setProperty("port", nPort.toString());

        Injector injector = Guice.createInjector(
                new ClothoTestModule(properties), 
                new MongoDBModule());
        
        SecurityManager securityManager = injector.getInstance(SecurityManager.class);
        SecurityUtils.setSecurityManager(securityManager);
        
        Persistor persistor = injector.getInstance(Persistor.class);
        persistor.deleteAll();
        TestUtils.setupTestData(persistor);
        
        ClothoWebserver server = injector.getInstance(ClothoWebserver.class);
        
        server.start();
    }
}

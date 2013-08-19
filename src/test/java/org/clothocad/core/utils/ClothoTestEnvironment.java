/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.utils;

import com.google.inject.Guice;
import com.google.inject.Injector;
import java.util.Properties;
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

        Injector injector = Guice.createInjector(new ClothoTestModule(properties), new MongoDBModule());
        
        Persistor persistor = injector.getInstance(Persistor.class);
        persistor.deleteAll();
        TestUtils.setupTestData(persistor);
        
        ClothoWebserver server = injector.getInstance(ClothoWebserver.class);
    }
}

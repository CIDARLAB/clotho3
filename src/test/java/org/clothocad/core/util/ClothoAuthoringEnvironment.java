/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import com.google.inject.Guice;
import com.google.inject.Injector;
import java.util.Arrays;
import java.util.Properties;
import org.apache.shiro.SecurityUtils;
import org.clothocad.core.ClothoStarter;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.persistence.mongodb.MongoDBModule;
import org.clothocad.core.security.ClothoRealm;
import static org.clothocad.core.util.ClothoTestEnvironment.main;
import org.clothocad.webserver.jetty.ClothoWebserver;

/**
 *
 * @author spaige
 */
public class ClothoAuthoringEnvironment extends ClothoStarter {

    public static void main(String[] args)
            throws Exception {

        Integer nPort = 8080;
        if (args.length > 1) {
            nPort = Integer.parseInt(args[0]);
        }

        Properties properties = new Properties();
        properties.setProperty("port", nPort.toString());

        Injector injector = Guice.createInjector(
                new ClothoAuthoringModule(properties),
                new MongoDBModule());

        org.apache.shiro.mgt.SecurityManager securityManager = injector.getInstance(org.apache.shiro.mgt.SecurityManager.class);
        SecurityUtils.setSecurityManager(securityManager);

        //test-specific setup

        Persistor persistor = injector.getInstance(Persistor.class);
        ClothoRealm realm = injector.getInstance(ClothoRealm.class);
        //persistor.deleteAll();
        TestUtils.importTestJSON(persistor);
        TestUtils.setupTestUsers(realm);

        server = injector.getInstance(ClothoWebserver.class);

        server.start();
    }

    @Override
    public void start() throws Exception {
        System.out.println("starting AuthoringEnvironment with arguments " + Arrays.toString(context.getArguments()));
        main(context.getArguments());
    }
}

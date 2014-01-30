/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import com.google.inject.Guice;
import com.google.inject.Injector;
import java.util.Arrays;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.apache.shiro.SecurityUtils;
import org.apache.shiro.mgt.SecurityManager;
import org.clothocad.core.ClothoStarter;
import static org.clothocad.core.ClothoStarter.main;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.persistence.jongo.JongoModule;
import org.clothocad.core.security.ClothoRealm;
import org.clothocad.core.testers.ClothoTestModule;
import org.clothocad.webserver.jetty.ClothoWebserver;

/**
 * Starts a clotho server instance configured for testing
 *
 * @author spaige
 */
public class ClothoTestEnvironment extends ClothoStarter {

    public static void main(String[] args)
            throws Exception {
        try {
            CommandLine cmd = parseArgs(args);

            if (cmd.hasOption("help")) {
                printHelp();
                return;
            }
            //TODO: if keystorepass option passed w/o arg, prompt for password 
            Injector injector = Guice.createInjector(new ClothoTestModule(commandToProperties(cmd)), new JongoModule());

            Persistor persistor = injector.getInstance(Persistor.class);

            SecurityManager securityManager = injector.getInstance(SecurityManager.class);
            SecurityUtils.setSecurityManager(securityManager);

            //test-specific setup
            ClothoRealm realm = injector.getInstance(ClothoRealm.class);
            persistor.deleteAll();
            realm.deleteAll();
            TestUtils.setupTestData(persistor);
            TestUtils.setupTestUsers(realm);

            server = injector.getInstance(ClothoWebserver.class);
            server.start();
        } catch (ParseException e) {
            //TODO: customise message to include default values
            System.out.println(e.getMessage());
            printHelp();
        }
    }

    @Override
    public void start() throws Exception {
        System.out.println("starting TestEnvironment with arguments " + Arrays.toString(context.getArguments()));
        main(context.getArguments());
    }
}

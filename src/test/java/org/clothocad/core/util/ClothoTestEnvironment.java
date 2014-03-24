/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import com.google.inject.Guice;
import com.google.inject.Injector;
import java.util.Arrays;
import java.util.Properties;
import org.apache.shiro.mgt.SecurityManager;
import org.apache.shiro.SecurityUtils;
import org.clothocad.core.AbstractClothoStarter;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.persistence.jongo.JongoModule;
import org.clothocad.core.security.ClothoRealm;
import org.clothocad.core.testers.ClothoTestModule;

/**
 * Starts a clotho server instance configured for testing
 *
 * @author spaige
 */
public class ClothoTestEnvironment extends AbstractClothoStarter {
    public static void main(String[] args) throws Exception {
        baseMain(args, new MainHook() {
            @Override public Injector
            getInjector(Properties config) {
                final Properties override = new Properties(config);
                override.setProperty("dbname", "testClotho");
                return Guice.createInjector(
                    new ClothoTestModule(override),
                    new JongoModule()
                );
            }

            @Override public void
            call(Injector injector) {
                SecurityManager securityManager =
                    injector.getInstance(SecurityManager.class);
                SecurityUtils.setSecurityManager(securityManager);

                //test-specific setup
                Persistor persistor = injector.getInstance(Persistor.class);
                ClothoRealm realm = injector.getInstance(ClothoRealm.class);
                persistor.deleteAll();
                realm.deleteAll();
                TestUtils.setupTestData(persistor);
                TestUtils.setupTestUsers(realm);
            }
        });
    }

    @Override
    public void start() throws Exception {
        System.out.println("starting TestEnvironment with arguments " + Arrays.toString(context.getArguments()));
        main(context.getArguments());
    }
}

package org.clothocad.core.util;

import com.google.inject.Guice;
import com.google.inject.Injector;
import java.util.Arrays;
import org.apache.commons.cli.CommandLine;
import org.apache.shiro.mgt.SecurityManager;
import org.apache.shiro.SecurityUtils;
import org.clothocad.core.AbstractClothoStarter;
import org.clothocad.core.persistence.mongodb.MongoDBModule;
import org.clothocad.core.persistence.Persistor;
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
            getInjector(CommandLine cmd) {
                return Guice.createInjector(
                    new ClothoTestModule(commandToProperties(cmd)),
                    new MongoDBModule()
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

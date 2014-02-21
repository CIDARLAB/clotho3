package org.clothocad.core.util;

import static org.clothocad.core.util.ClothoTestEnvironment.main;

import com.google.inject.Guice;
import com.google.inject.Injector;
import com.google.inject.Key;
import com.google.inject.name.Names;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import org.apache.commons.cli.CommandLine;
import org.apache.shiro.mgt.SecurityManager;
import org.apache.shiro.SecurityUtils;
import org.clothocad.core.AbstractClothoStarter;
import org.clothocad.core.persistence.dataauthoring.FileHookPersistor;
import org.clothocad.core.persistence.mongodb.MongoDBModule;
import org.clothocad.core.security.ClothoRealm;

/**
 *
 * @author spaige
 */
public class ClothoAuthoringEnvironment extends AbstractClothoStarter {
    public static void main(String[] args) throws Exception {
        baseMain(args, new MainHook() {
            @Override public Injector
            getInjector(CommandLine cmd) {
                return Guice.createInjector(
                    new ClothoAuthoringModule(commandToProperties(cmd)),
                    new MongoDBModule());
            }

            @Override public void
            call(Injector injector) {
                SecurityManager securityManager
                    = injector.getInstance(SecurityManager.class);
                SecurityUtils.setSecurityManager(securityManager);

                /* test-specific setup */
                FileHookPersistor persistor =
                    injector.getInstance(FileHookPersistor.class);
                ClothoRealm realm = injector.getInstance(ClothoRealm.class);
                Path storageFolder = injector.getInstance(
                    Key.get(Path.class, Names.named("storagefolder"))
                );
                if (isBoringDirectory(storageFolder))
                    persistor.initializeBuiltInSchemas();

                TestUtils.importTestJSON(
                    Paths.get("src", "test", "resources").toString(),
                    persistor,
                    false
                );
                TestUtils.importTestJSON(
                    storageFolder.toString(), persistor, true);
                TestUtils.importTestJSON(persistor);
                TestUtils.setupTestUsers(realm);
            }

            private boolean isBoringDirectory(Path dir) {
                try {
                    return !Files.newDirectoryStream(dir).iterator().hasNext();
                } catch (IOException e) {
                    return true;
                }
            }
        });
    }

    @Override
    public void start() throws Exception {
        System.out.println("starting AuthoringEnvironment with arguments " + Arrays.toString(context.getArguments()));
        main(context.getArguments());
    }
}

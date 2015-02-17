package org.clothocad.core.util;

import static org.clothocad.core.util.ClothoTestEnvironment.main;

import org.clothocad.core.AbstractClothoStarter;
import org.clothocad.core.persistence.dataauthoring.FileHookPersistor;
import org.clothocad.core.persistence.jongo.JongoModule;
import org.clothocad.core.security.ClothoRealm;

import com.google.inject.Guice;
import com.google.inject.Injector;
import com.google.inject.Key;
import com.google.inject.name.Names;

import org.apache.shiro.mgt.SecurityManager;
import org.apache.shiro.SecurityUtils;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Properties;

/**
 *
 * @author spaige
 * @author billcao
 */
public class ClothoAuthoringEnvironment extends AbstractClothoStarter {
    public static void main(String[] args) throws Exception {
        baseMain(args, new MainHook() {
            @Override public Injector
            getInjector(Properties config) {
                final Properties override = new Properties(config);
                override.setProperty("dbname", "authoringenv");
                return Guice.createInjector(
                    new ClothoAuthoringModule(override),
                    new JongoModule()
                );
            }

            @Override public void
            call(Injector injector) {
                SecurityManager securityManager
                    = injector.getInstance(SecurityManager.class);
                SecurityUtils.setSecurityManager(securityManager);

                /* test-specific setup */
                FileHookPersistor persistor =
                    injector.getInstance(FileHookPersistor.class);

                Path storageFolder = injector.getInstance(
                    Key.get(Path.class, Names.named("storagefolder"))
                );
                Path resourceFolder = injector.getInstance(
                    Key.get(Path.class, Names.named("resourcefolder"))
                );

                persistor.initializeBuiltInSchemas();

                TestUtils.importTestJSON(
                    storageFolder, persistor.getConnection(), true);
                TestUtils.importResources(
                    resourceFolder, persistor.getConnection(), true);
            }
        });
    }

    @Override
    public void start() throws Exception {
        System.out.println("starting AuthoringEnvironment with arguments " + Arrays.toString(context.getArguments()));
        main(context.getArguments());
    }
}

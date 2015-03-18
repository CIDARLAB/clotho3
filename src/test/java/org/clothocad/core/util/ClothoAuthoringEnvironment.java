package org.clothocad.core.util;

import org.clothocad.core.security.nosecurity.NoSecurityModule;
import static org.clothocad.core.util.ClothoTestEnvironment.main;

import com.google.inject.Guice;
import com.google.inject.Injector;
import com.google.inject.Key;
import com.google.inject.name.Names;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Properties;
import javax.servlet.ServletRequest;
import javax.servlet.ServletResponse;
import org.apache.shiro.SecurityUtils;
import org.apache.shiro.mgt.DefaultSecurityManager;
import org.apache.shiro.mgt.SecurityManager;
import org.apache.shiro.session.Session;
import org.apache.shiro.subject.PrincipalCollection;
import org.apache.shiro.subject.Subject;
import org.apache.shiro.subject.SubjectContext;
import org.apache.shiro.web.mgt.DefaultWebSubjectFactory;
import org.apache.shiro.web.subject.WebSubjectContext;
import org.clothocad.core.AbstractClothoStarter;
import org.clothocad.core.persistence.dataauthoring.FileHookPersistor;
import org.clothocad.core.persistence.jongo.JongoModule;
import org.clothocad.core.security.ClothoRealm;

/**
 *
 * @author spaige
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
                    new NoSecurityModule(),
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
                ClothoRealm realm = injector.getInstance(ClothoRealm.class);
                Path storageFolder = injector.getInstance(
                    Key.get(Path.class, Names.named("storagefolder"))
                );
                
                //clear db so that authored JSON directory is single source of truth
                persistor.deleteAll();
                persistor.initializeBuiltInSchemas();

                TestUtils.importJSONFromDirectory(
                    storageFolder, persistor.getConnection(), null, true, false);
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

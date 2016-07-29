package org.clothocad.core.util;

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
                    new AuthoringSecurityModule(),
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

//                TestUtils.importJSONFromDirectory(
//                    storageFolder, persistor.getConnection(), null, true, false);
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
    
    public static class AuthoringSubjectFactory extends DefaultWebSubjectFactory {

        @Override
        public Subject createSubject(SubjectContext context) {
            if (!(context instanceof WebSubjectContext)) {
                SecurityManager securityManager = context.resolveSecurityManager();
                Session session = context.resolveSession();
                boolean sessionCreationEnabled = context.isSessionCreationEnabled();
                PrincipalCollection principals = context.resolvePrincipals();
                boolean authenticated = context.resolveAuthenticated();
                String host = context.resolveHost();

                return new LoggedInSubject(principals, authenticated, host, session, sessionCreationEnabled, securityManager);
            }
            WebSubjectContext wsc = (WebSubjectContext) context;
            SecurityManager securityManager = wsc.resolveSecurityManager();
            Session session = wsc.resolveSession();
            boolean sessionEnabled = wsc.isSessionCreationEnabled();
            PrincipalCollection principals = wsc.resolvePrincipals();
            boolean authenticated = wsc.resolveAuthenticated();
            String host = wsc.resolveHost();
            ServletRequest request = wsc.resolveServletRequest();
            ServletResponse response = wsc.resolveServletResponse();

            return new LoggedInWebSubject(principals, authenticated, host, session, sessionEnabled,
                    request, response, securityManager);
        }
    }
}

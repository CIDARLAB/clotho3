/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import com.google.inject.binder.AnnotatedBindingBuilder;
import javax.inject.Inject;
import javax.inject.Named;
import org.apache.shiro.guice.web.ShiroWebModule;
import org.apache.shiro.session.mgt.SessionManager;
import org.apache.shiro.web.session.mgt.DefaultWebSessionManager;
import org.eclipse.jetty.servlet.ServletContextHandler;

/**
 *
 * @author spaige
 */
public class SecurityModule extends ShiroWebModule {

    //XXX: augh jetty hack
    @Inject
    public SecurityModule(@Named("containerServletContext") ServletContextHandler servletContext) {
        super(servletContext.getServletContext());
    }


    @Override
    protected void configureShiroWeb() {
            bindRealm().to(ClothoRealm.class);
            
            ShiroWebModule.bindGuiceFilter(binder());
    }

    @Override
    protected void bindSessionManager(AnnotatedBindingBuilder<SessionManager> bind) {
        bind.to(DefaultWebSessionManager.class);
    }
    
    
}

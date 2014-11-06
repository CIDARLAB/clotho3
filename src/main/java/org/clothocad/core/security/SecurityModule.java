/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import com.google.inject.binder.AnnotatedBindingBuilder;
import com.google.inject.name.Names;
import org.apache.shiro.guice.web.ShiroWebModule;
import org.apache.shiro.session.mgt.SessionManager;
import org.apache.shiro.web.session.mgt.DefaultWebSessionManager;
import org.eclipse.jetty.servlet.ServletContextHandler;

/**
 *
 * @author spaige
 */
public class SecurityModule extends ShiroWebModule {

    protected ServletContextHandler handler;
    
    public SecurityModule(ServletContextHandler servletContextHandler) {
        super(servletContextHandler.getServletContext());
        handler = servletContextHandler;
    }

    public SecurityModule(){
        this(new ServletContextHandler());
    }

    @Override
    protected void configureShiroWeb() {
            bindRealm().to(ClothoRealm.class);
            //bindRealm().to(OAuthRealm.class);
            bind(ServletContextHandler.class).annotatedWith(Names.named("containerServletContext")).toInstance(handler);
            expose(ServletContextHandler.class).annotatedWith(Names.named("containerServletContext"));
            ShiroWebModule.bindGuiceFilter(binder());
    }

    @Override
    protected void bindSessionManager(AnnotatedBindingBuilder<SessionManager> bind) {
        bind.to(DefaultWebSessionManager.class);
    }
}

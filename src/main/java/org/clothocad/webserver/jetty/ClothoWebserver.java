package org.clothocad.webserver.jetty;

import org.clothocad.core.communication.Router;
import org.clothocad.core.communication.ws.ClothoWebSocket;

import com.google.inject.servlet.GuiceFilter;

import lombok.Getter;

import org.eclipse.jetty.security.ConstraintMapping;
import org.eclipse.jetty.security.ConstraintSecurityHandler;
import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.server.handler.HandlerList;
//import org.eclipse.jetty.server.nio.SelectChannelConnector;
//import org.eclipse.jetty.server.ssl.SslConnector;
import org.eclipse.jetty.servlet.DefaultServlet;
import org.eclipse.jetty.servlet.ServletContextHandler;
import org.eclipse.jetty.servlet.ServletHolder;
import org.eclipse.jetty.util.security.Constraint;
import org.eclipse.jetty.websocket.servlet.WebSocketServlet;
import org.eclipse.jetty.websocket.servlet.WebSocketServletFactory;

import javax.inject.Inject;
import javax.inject.Named;
import javax.servlet.http.HttpServletRequest;
import org.eclipse.jetty.server.HttpConfiguration;
import org.eclipse.jetty.server.HttpConnectionFactory;
import org.eclipse.jetty.server.SecureRequestCustomizer;
import org.eclipse.jetty.server.ServerConnector;
import org.eclipse.jetty.util.ssl.SslContextFactory;

//TODO: convert config to guice module
//TODO: make easy to switch ssl requirement on/off for deploy testing
public class ClothoWebserver {

    @Getter
    final Server server;

    @Inject
    public ClothoWebserver(@Named("port") int nPort,
            @Named("confidentialport") int confidentialPort,
            SslContextFactory sslContextFactory,
            @Named("containerServletContext") ServletContextHandler servletHandler,
            final Router router, @Named("clientdirectory") String clientDirectory)
            throws Exception {

        server = new Server();

        //Connectors
        
        HttpConfiguration http_config = new HttpConfiguration();
        http_config.setSecurePort(confidentialPort);
        http_config.setOutputBufferSize(8192);
        
        
        ServerConnector http = new ServerConnector(server, new HttpConnectionFactory(http_config));
        http.setPort(nPort);
        http.setIdleTimeout(3600000);
        
        HttpConfiguration https_config = new HttpConfiguration(http_config);
        SecureRequestCustomizer src = new SecureRequestCustomizer();
        src.setStsMaxAge(2000);
        src.setStsIncludeSubDomains(true);
        https_config.setOutputBufferSize(8192);
        https_config.addCustomizer(src);
        
        ServerConnector https = new ServerConnector(server, sslContextFactory, new HttpConnectionFactory(https_config));
        https.setPort(confidentialPort);
        https.setIdleTimeout(3600000);
        
        server.addConnector(http);
        server.addConnector(https);

        // Connection constraints
        Constraint constraint = new Constraint();
        constraint.setDataConstraint(Constraint.DC_CONFIDENTIAL);

        ConstraintMapping cm = new ConstraintMapping();
        cm.setConstraint(constraint);
        cm.setPathSpec("/*");

        ConstraintSecurityHandler constraintHandler = new ConstraintSecurityHandler();
        constraintHandler.setConstraintMappings(new ConstraintMapping[]{cm});


        // Static resources
        DefaultServlet staticServlet = new DefaultServlet();

        // Handler stack
        servletHandler.setContextPath("/");
        servletHandler.setResourceBase(clientDirectory);
        servletHandler.setWelcomeFiles(new String[]{"index.html"});

        servletHandler.addFilter(GuiceFilter.class, "/*", null);

        servletHandler.addServlet(new ServletHolder(staticServlet), "/*");
        servletHandler.addServlet(new ServletHolder(ClothoServlet.class), "/websocket");
        servletHandler.addServlet(new ServletHolder(new RestApi(router)), "/data/*");

        HandlerList handlers = new HandlerList();
        handlers.addHandler(constraintHandler);
        constraintHandler.setHandler(servletHandler);
        server.setHandler(handlers);
    }
    
    @SuppressWarnings("serial")
    public class ClothoServlet extends WebSocketServlet
    {
        @Override
        public void configure(WebSocketServletFactory factory)
        {
            factory.register(ClothoWebSocket.class);
        }
    }

    public void start() throws Exception {
        
        server.start();
        //server.doStart();
    }
}

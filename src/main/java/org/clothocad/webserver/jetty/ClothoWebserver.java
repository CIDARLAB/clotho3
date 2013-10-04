package org.clothocad.webserver.jetty;

import com.google.inject.servlet.GuiceFilter;
import java.security.KeyStore;
import javax.inject.Inject;
import javax.inject.Named;
import javax.servlet.http.HttpServletRequest;
import lombok.Getter;
import org.clothocad.core.layers.communication.Router;

import org.clothocad.core.layers.communication.connection.ws.ClothoWebSocket;
import org.eclipse.jetty.security.ConstraintMapping;
import org.eclipse.jetty.security.ConstraintSecurityHandler;
import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.server.handler.HandlerList;
import org.eclipse.jetty.server.nio.SelectChannelConnector;
import org.eclipse.jetty.websocket.WebSocket;
import org.eclipse.jetty.server.ssl.SslSelectChannelConnector;
import org.eclipse.jetty.servlet.DefaultServlet;
import org.eclipse.jetty.servlet.ServletContextHandler;
import org.eclipse.jetty.servlet.ServletHolder;
import org.eclipse.jetty.util.security.Constraint;
import org.eclipse.jetty.util.ssl.SslContextFactory;
import org.eclipse.jetty.websocket.WebSocketServlet;

//TODO: convert config to guice module
public class ClothoWebserver {

    @Getter
    final Server server;
    

    
    @Inject
    public ClothoWebserver(@Named("port") int nPort,
            KeyStore keystore, 
            @Named("containerServletContext") ServletContextHandler servletHandler,
            final Router router)
            throws Exception {

        int confidentialPort = 8443; //TODO: make configurable
        
        server = new Server();

        //Connectors

        SelectChannelConnector connector0 = new SelectChannelConnector();
        connector0.setPort(nPort);
        connector0.setMaxIdleTime(3600000);
        connector0.setRequestHeaderSize(8192);
        connector0.setConfidentialPort(confidentialPort);
        server.addConnector(connector0);
        
        SslSelectChannelConnector ssl_connector = new SslSelectChannelConnector();
        ssl_connector.setPort(confidentialPort); 
        ssl_connector.setMaxIdleTime(3600000);
        SslContextFactory cf = ssl_connector.getSslContextFactory();
        cf.setKeyStore(keystore);
        server.addConnector(ssl_connector);
        
        
        //Connection constraints
        
        Constraint constraint = new Constraint();
        constraint.setDataConstraint(Constraint.DC_CONFIDENTIAL);

        ConstraintMapping cm = new ConstraintMapping();
        cm.setConstraint(constraint);
        cm.setPathSpec("/*");
        
        ConstraintSecurityHandler constraintHandler = new ConstraintSecurityHandler();
        constraintHandler.setConstraintMappings(new ConstraintMapping[]{cm});
        
        //Websocket

        WebSocketServlet wsServlet = new WebSocketServlet() {

            @Override
            public WebSocket doWebSocketConnect(HttpServletRequest request, String protocol) {
                return new ClothoWebSocket(request.getSession().getId(), router); 
            }
            
        };
                
        //Static resources
        
        DefaultServlet staticServlet = new DefaultServlet();
        
        //Handler stack
        
        
        
        servletHandler.setContextPath("/");
        servletHandler.setResourceBase("./clotho3-web/"); 
        servletHandler.setWelcomeFiles(new String[]{"index.html"});
        
        servletHandler.addFilter(GuiceFilter.class, "/*", null);
        servletHandler.addServlet(new ServletHolder(staticServlet), "/*");
        servletHandler.addServlet(new ServletHolder(wsServlet), "/websocket");
                
        HandlerList handlers = new HandlerList();
        handlers.addHandler(constraintHandler);
        constraintHandler.setHandler(servletHandler);        
        server.setHandler(handlers);

    }
    
    
    public void start() throws Exception {
        server.start();
    }
}
package org.clothocad.webserver.jetty;

import com.google.inject.servlet.GuiceFilter;
import java.security.KeyStore;
import javax.inject.Inject;
import javax.inject.Named;
import javax.servlet.http.HttpServletRequest;
import lombok.Getter;

import org.clothocad.core.layers.communication.connection.ws.ClothoWebSocket;
import org.eclipse.jetty.server.Connector;
import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.server.nio.SelectChannelConnector;
import org.eclipse.jetty.websocket.WebSocket;
import org.eclipse.jetty.server.ssl.SslSelectChannelConnector;
import org.eclipse.jetty.servlet.DefaultServlet;
import org.eclipse.jetty.servlet.ServletContextHandler;
import org.eclipse.jetty.servlet.ServletHolder;
import org.eclipse.jetty.util.ssl.SslContextFactory;
import org.eclipse.jetty.websocket.WebSocketServlet;

public class ClothoWebserver {

    @Getter
    final Server server;

    @Inject
    public ClothoWebserver(@Named("port") int nPort,
            KeyStore keystore)
            throws Exception {

        server = new Server(nPort);

        //Connectors

        SelectChannelConnector connector0 = new SelectChannelConnector();
        connector0.setPort(nPort);
        connector0.setMaxIdleTime(30000);
        connector0.setRequestHeaderSize(8192);

        SslSelectChannelConnector ssl_connector = new SslSelectChannelConnector();
        ssl_connector.setPort(8443); //TODO: make configurable
        SslContextFactory cf = ssl_connector.getSslContextFactory();
        cf.setKeyStore(keystore);
        server.setConnectors(new Connector[]{connector0, ssl_connector});
        
        //Websocket

        WebSocketServlet wsServlet = new WebSocketServlet() {

            @Override
            public WebSocket doWebSocketConnect(HttpServletRequest request, String protocol) {
                return new ClothoWebSocket(request.getSession().getId()); //todo - can inject router here w/ guice
            }
            
        };
        
        //Static resources
        
        DefaultServlet staticServlet = new DefaultServlet();
        
        //Handler stack
        
        ServletContextHandler servletHandler = new ServletContextHandler(ServletContextHandler.SESSIONS);
        servletHandler.setContextPath("/");
        servletHandler.setResourceBase("./clotho3-web/");
        servletHandler.setWelcomeFiles(new String[]{"index.html"});
        
        servletHandler.addFilter(GuiceFilter.class, "/*", null);
        servletHandler.addServlet(new ServletHolder(staticServlet), "/*");
        servletHandler.addServlet(new ServletHolder(wsServlet), "/websocket");
        
        server.setHandler(servletHandler);

        server.start();
    }
}
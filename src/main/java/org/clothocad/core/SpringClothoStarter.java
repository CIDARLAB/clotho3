/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core;

import org.clothocad.core.communication.Router;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.security.ClothoRealm;
import org.clothocad.webserver.jetty.ClothoWSHandler;
import org.eclipse.jetty.websocket.api.WebSocketBehavior;
import org.eclipse.jetty.websocket.api.WebSocketPolicy;
import org.eclipse.jetty.websocket.server.WebSocketServerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.boot.builder.SpringApplicationBuilder;
import org.springframework.boot.context.web.SpringBootServletInitializer;
import org.springframework.boot.SpringApplication;
import org.springframework.context.ConfigurableApplicationContext;
import org.springframework.context.annotation.Bean;
import org.springframework.web.socket.config.annotation.EnableWebSocket;
import org.springframework.web.socket.config.annotation.WebSocketConfigurer;
import org.springframework.web.socket.config.annotation.WebSocketHandlerRegistry;
import org.springframework.web.socket.server.HandshakeHandler;
import org.springframework.web.socket.server.jetty.JettyRequestUpgradeStrategy;
import org.springframework.web.socket.server.support.DefaultHandshakeHandler;

@EnableWebSocket
@SpringBootApplication
public class SpringClothoStarter extends SpringBootServletInitializer implements WebSocketConfigurer {

    static ConfigurableApplicationContext ctx;
    static Router rbean;
    static Persistor pbean;
    static ClothoRealm crbean;
    static ClothoWSHandler clothoPls;
    
    @Override
    protected SpringApplicationBuilder configure(SpringApplicationBuilder application) {
        return application.sources(SpringClothoStarter.class);
    }

    @Override
    public void registerWebSocketHandlers(WebSocketHandlerRegistry registry) {        
        clothoPls = new ClothoWSHandler();        
        registry.addHandler(clothoPls, "/websocket").setHandshakeHandler(factoryBean());
    }
    
    //Configure buffer size and timeouts
    @Bean
    public HandshakeHandler factoryBean()
    {
        WebSocketPolicy policy = new WebSocketPolicy(WebSocketBehavior.SERVER);
        policy.setInputBufferSize(8192);
        policy.setIdleTimeout(3600000);
        WebSocketServerFactory factory = new WebSocketServerFactory(policy);
        return new DefaultHandshakeHandler(new JettyRequestUpgradeStrategy(factory));
    }
    

    public static void main(String[] args) throws Exception {
        ctx = SpringApplication.run(SpringClothoStarter.class, args);
        
        rbean = ctx.getBean("router", Router.class);
        pbean = ctx.getBean("persistor", Persistor.class);
        
        clothoPls.setRouter(rbean);
        
        System.out.println("rbean: " + rbean.toString());
        System.out.println("pbean: " + pbean.toString());
    }
}

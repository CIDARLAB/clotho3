/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.webserver.jetty;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.context.ApplicationListener;
import org.springframework.stereotype.Component;

/**
 *
 * @author david
 */
@Component
public class ApplicationEnvironmentPreparedEvent implements ApplicationListener<org.springframework.boot.context.event.ApplicationEnvironmentPreparedEvent> {

    public RouterBean router;
    
    @Autowired
    public void setRouterBean(RouterBean router)
    {
        this.router = router;
    }
    
    @Override
    public void onApplicationEvent(org.springframework.boot.context.event.ApplicationEnvironmentPreparedEvent event) {
        System.out.println();
        System.out.println();
        System.out.println();
        System.out.println("At the ApplicationEnvironmentPreparedEvent for the router bean, hello hello I'm here");
        System.out.println();
        System.out.println();
        System.out.println();
        
        router.init();
    }
    
}

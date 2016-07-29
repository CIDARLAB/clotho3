/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.webserver.jetty;

import org.clothocad.core.ClothoModule;
import org.clothocad.core.communication.Router;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.persistence.jongo.JongoModule;
import org.clothocad.core.security.SecurityModule;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

/**
 *
 * @author david
 */
@Component
public class RouterBean {
    
    public ClothoModule clothoMod;
    public SecurityModule securityMod;
    public JongoModule jongoMod;
    public boolean yes = false;
    
    @Autowired
    public Persistor persistor;
    
    public void init()
    {
        clothoMod = new ClothoModule(null);
        securityMod = new SecurityModule();
        jongoMod = new JongoModule();
        yes = true;
    }
    
    
}

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import org.clothocad.core.communication.Router;
import org.clothocad.core.communication.ServerSideAPI;
import org.clothocad.core.communication.TestConnection;
import org.clothocad.core.execution.Mind;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.util.SecurityTestUtils;
import org.junit.After;
import org.junit.Before;

/**
 *
 * @author spaige
 */
public abstract class AbstractSecurityTest extends AbstractShiroTest {
    protected Persistor persistor;
    protected ServerSideAPI api;
    protected SecurityTestUtils util;

    public AbstractSecurityTest() {
        super();
        persistor = injector.getInstance(Persistor.class);
        persistor.deleteAll();
        util = new SecurityTestUtils(persistor, injector.getInstance(ClothoRealm.class));
    
    }
    
    @Before
    @Override
    public void beforeTest() {
        super.beforeTest(); //To change body of generated methods, choose Tools | Templates.
    }

    @After
    @Override
    public void afterTest() {
        if (api != null) api.logout();
        super.afterTest(); //To change body of generated methods, choose Tools | Templates.
    }    
    
    
    /**
     * create a new instance of ServerSideAPI
     *
     * @param id String format id of ServerSideAPI
     */
    public void initAPI(String id) {
        Mind mind = new Mind();
        mind.setConnection(new TestConnection(id));
        api = new ServerSideAPI(mind, persistor, injector.getInstance(Router.class), id);
    }    

    public void login(String account){
        api.login(account, account);
    }
}



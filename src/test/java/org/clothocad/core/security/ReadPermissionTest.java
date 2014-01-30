/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import com.google.inject.Guice;
import com.google.inject.Injector;
import java.util.HashMap;
import java.util.Map;
import org.apache.shiro.SecurityUtils;
import org.apache.shiro.UnavailableSecurityManagerException;
import org.apache.shiro.authc.*;
import org.apache.shiro.authz.*;

import org.apache.shiro.subject.Subject;
import org.clothocad.core.communication.Router;
import org.clothocad.core.communication.ServerSideAPI;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.persistence.jongo.JongoModule;
import org.clothocad.core.testers.ClothoTestModule;
import org.clothocad.core.util.SecurityTestUtils;
import org.junit.Test;

/**
 * test case of user with 'read' permission test actions of read, edit, delete
 * and edit permission
 *
 * @author yu
 * @version 1.0
 */
public class ReadPermissionTest {

    private Router router;
    private Persistor persistor;
    private Injector injector;
    private ServerSideAPI api;
    private SecurityTestUtils util;

    /**
     * constructor
     */
    public ReadPermissionTest() {
    }

    /**
     * create a new instance of ServerSideAPI
     *
     * @param id String format id of ServerSideAPI
     */
    public void initAPI(String id) {
        injector = Guice.createInjector(new ClothoTestModule(), new JongoModule());
        persistor = injector.getInstance(Persistor.class);
        //ServerSideAPI api = new DummyAPI(persistor);
        api = new ServerSideAPI(null, persistor, null, id);
        util = new SecurityTestUtils(persistor);
    }

    /**
     * create a private object in ServerSideAPI
     *
     * @return returns the default object in ServerSideAPI
     */
    public Map<String, Object> initObj() {
        Map<String, Object> defObj = new HashMap<>();
        Subject currentUser = SecurityUtils.getSubject();
        //set object as private
        defObj.put("private", util.objects.get("private"));
        UsernamePasswordToken token = new UsernamePasswordToken("owner", "owner");
        currentUser.login(token);
        //log into ServerSideAPI to create an object
        api.login("owner", "owner");
        api.create(defObj);
        api.logout();
        return defObj;
    }

    /**
     * test read permission
     * @exception no exception expected
     */
    @Test
    public void testRead() {
        initAPI("0000");
        Map<String, Object> newObj = initObj();

        try {
            Subject currentUser = SecurityUtils.getSubject();

            UsernamePasswordToken token = new UsernamePasswordToken("read", "read");
            currentUser.login(token);
            token.setRememberMe(true);
            api.login("read", "read");
            api.get(newObj);
        } catch (UnavailableSecurityManagerException e) {
        }

    }

    /**
     * test edit permission, expecting UnauthorizedException
     *
     * @exception UnauthorizedException expected
     */
    @Test(expected = UnauthorizedException.class)
    public void testEdit() {
        initAPI("0001");
        Map<String, Object> newObj = initObj();

        try {
            Subject currentUser = SecurityUtils.getSubject();
            UsernamePasswordToken token = new UsernamePasswordToken("read", "read");
            currentUser.login(token);
            token.setRememberMe(true);
            api.login("read", "read");
            api.set(newObj);
        } catch (UnavailableSecurityManagerException e) {
        }

    }

    /**
     * test delete permission, expecting UnauthorizedException
     *
     * @exception UnauthorizedException expected
     */
    @Test(expected = UnauthorizedException.class)
    public void testDelete() {
        initAPI("0002");
        Map<String, Object> newObj = initObj();

        try {
            //new user
            Subject currentUser = SecurityUtils.getSubject();
            UsernamePasswordToken token = new UsernamePasswordToken("read", "read");
            currentUser.login(token);
            token.setRememberMe(true);
            //log in
            api.login("read", "read");
            //delete
            api.destroy(newObj);
        } catch (UnavailableSecurityManagerException e) {
        }

    }

    /**
     * test editing other user's permission, expecting UnauthorizedException
     *
     * @exception UnauthorizedException expected
     */
    @Test(expected = UnauthorizedException.class)
    public void testEditPermission() {
        initAPI("0003");
        Map<String, Object> newObj = initObj();

        try {
            Subject currentUser = SecurityUtils.getSubject();
            UsernamePasswordToken token = new UsernamePasswordToken("unauth", "unauth");
            currentUser.login(token);
            token.setRememberMe(true);
            api.login("unauth", "unauth");
            /*
             *code here to edit permission 
             */

        } catch (UnavailableSecurityManagerException e) {
        }

    }
}
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
import org.clothocad.core.persistence.mongodb.MongoDBModule;
import org.clothocad.core.testers.ClothoTestModule;
import org.clothocad.core.util.SecurityTestUtils;

import org.junit.Test;

/**
 * test case of user with 'none' permission test actions of read, edit, delete
 * and edit permission on a private object
 *
 * @author yu
 * @version 1.0
 */
public class UnauthPrivateTest {

    private Router router;
    private Persistor persistor;
    private Injector injector;
    private ServerSideAPI api;
    private SecurityTestUtils util;

    /**
     * constructor
     */
    public UnauthPrivateTest() {
    }

    /**
     * create a new instance of ServerSideAPI
     *
     * @param id String format id of ServerSideAPI
     */
    public void initAPI(String id) {
        injector = Guice.createInjector(new ClothoTestModule(), new MongoDBModule());
        persistor = injector.getInstance(Persistor.class);
        //ServerSideAPI api = new DummyAPI(persistor);
        api = new ServerSideAPI(null, persistor, null, id);
    }

    /**
     * create a new instance of ServerSideAPI
     *
     * @return returns the default object in ServerSideAPI
     */
    public Map<String, Object> initObj() {
        Map<String, Object> defObj = new HashMap<>();
        Subject currentUser = SecurityUtils.getSubject();
        currentUser.checkPermission(util.credentials.get("owner"));
        defObj.put("private", util.objects.get("private"));

        UsernamePasswordToken token = new UsernamePasswordToken("owner", "owner");
        currentUser.login(token);
        api.login("owner", "owner");
        api.create(defObj);
        api.logout();
        return defObj;
    }

    /**
     * test read permission, should catch UnauthorizedException
     *
     * @exception UnauthorizedException expected
     */
    @Test(expected = UnauthorizedException.class)
    public void testRead() {
        initAPI("0000");
        Map<String, Object> newObj = initObj();

        try {
            //SecurityManager securityManager = SecurityUtils.getSecurityManager();
            //Subject currentUser = new Subject.Builder(securityManager).buildSubject();
            Subject currentUser = SecurityUtils.getSubject();
            UsernamePasswordToken token = new UsernamePasswordToken("none", "none");
            currentUser.login(token);
            token.setRememberMe(true);
            api.login("none", "none");
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
            UsernamePasswordToken token = new UsernamePasswordToken("none", "none");
            currentUser.login(token);
            token.setRememberMe(true);
            api.login("none", "none");
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
            Subject currentUser = SecurityUtils.getSubject();
            UsernamePasswordToken token = new UsernamePasswordToken("none", "none");
            currentUser.login(token);
            token.setRememberMe(true);
            api.login("none", "none");
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
            UsernamePasswordToken token = new UsernamePasswordToken("none", "none");
            currentUser.login(token);
            token.setRememberMe(true);
            api.login("none", "none");
            /*
             *code here to edit permission 
             */

        } catch (UnavailableSecurityManagerException e) {
        }

    }
}
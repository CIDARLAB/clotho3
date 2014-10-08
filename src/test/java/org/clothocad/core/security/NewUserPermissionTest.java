/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import javax.persistence.EntityNotFoundException;
import javax.script.ScriptException;
import org.apache.shiro.authz.*;

import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import static org.clothocad.core.security.ClothoPermission.WRITE;
import org.clothocad.model.Institution;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 * test case of user with 'write' permission test actions of read, edit, delete
 * and edit permission on a private object
 *
 * @author yu
 * @version 1.0
 */
public class NewUserPermissionTest extends AbstractSecurityTest {

    /**
     * constructor
     */
    public NewUserPermissionTest() {
        super();
    }

    /**
     * test read action
     *
     */
    @Test
    public void testReadPublic() {
        initAPI("0000");

         ObjBase pub = util.getPublic();
        api.login("none", "none");
        assertEquals(pub.getId(), persistor.get(Institution.class, pub.getId()).getId());
    }

    @Test(expected = EntityNotFoundException.class)
    public void testReadPrivate() {
        initAPI("0000");

        ObjBase priv = util.getPrivate();
        api.login("none", "none");
        assertEquals(priv.getId(), persistor.get(Institution.class, priv.getId()).getId());
    }

    @Test
    public void testRunPublic() throws ScriptException, IllegalAccessException, IllegalArgumentException, InvocationTargetException{
        initAPI("runRun");
        login("none");
        Map<String,Object> command = new HashMap<>();
        command.put("args", new ArrayList());
        command.put("id", util.getPublicModule().getId());
        command.put("function", "function");
        assertEquals("function ran!", api.run(command));
    }

    @Test(expected = EntityNotFoundException.class)
    public void testRunPrivate() throws ScriptException, IllegalAccessException, IllegalArgumentException, InvocationTargetException{
        initAPI("runRun");
        login("none");
        Map<String,Object> command = new HashMap<>();
        command.put("args", new ArrayList());
        command.put("id", util.getPrivateModule().getId());
        command.put("function", "function");
        assertEquals("function ran!", api.run(command));
        }


    /**
     * test edit action
     *
     */
    @Test(expected = AuthorizationException.class)
    public void testEdit() {
        initAPI("0001");
        ObjBase priv = util.getPrivate();
        try{

            Map<String,Object> edit = new HashMap<>();
            edit.put("id", priv.getId());
            edit.put("name", "Changed Name");
            api.login("none", "none");
            persistor.save(edit);
        } finally {
            //reset state of persistor
            persistor.save(priv);
        }

    }

    /**
     * test delete action
     *
     * @exception UnauthorizedException
     */
    @Test(expected = UnauthorizedException.class)
    public void testDelete() {
        initAPI("0002");
            ObjectId id = util.getPrivate().getId();
            api.login("none", "none");
            persistor.delete(id);
    }
    /**
     * test edit permission
     *
     * @exception UnauthorizedException
     */
    @Test(expected = AuthorizationException.class)
    public void testGrant() {
        initAPI("0003");

        api.login("none", "none");
        //private
        realm.addPermissions(ClothoRealm.ANONYMOUS_USER, WRITE.actions, util.getPrivate().getId());
        //public
        realm.addPermissions(ClothoRealm.ANONYMOUS_USER, WRITE.actions, util.getPrivate().getId());        
    }
}
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import com.google.common.collect.Lists;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
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
public class WritePermissionTest extends AbstractSecurityTest {

    /**
     * constructor
     */
    public WritePermissionTest() {
        super();
    }

    /**
     * test read action
     *
     */
    @Test
    public void testRead() {
        initAPI("0000");

        ObjBase priv = util.getPrivate();
            api.login("write", "write");
        assertEquals(priv.getId(), persistor.get(Institution.class, priv.getId()).getId());
    }
    
    @Test
    public void testFind() {
        initAPI("writeFind");
        api.login("write", "write");
        Map<String, Object> query = new HashMap<>();
        query.put("schema", Institution.class.getName());
        assertEquals(2, Lists.newArrayList(persistor.find(query)).size());
        }

    @Test
    public void testRun() throws Exception {
        initAPI("writeRun");
        login("owner");
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
    @Test
    public void testEdit() {
        initAPI("0001");
        ObjBase priv = util.getPrivate();
        try{

            Map<String,Object> edit = new HashMap<>();
            edit.put("id", priv.getId());
            edit.put("name", "Changed Name");
            api.login("write", "write");
            persistor.save(edit);
            assertEquals("Changed Name", persistor.get(Institution.class, priv.getId()).getName());
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
            api.login("write", "write");
            persistor.delete(id);
    }
    /**
     * test edit permission
     *
     * @exception UnauthorizedException
     */
    @Test(expected = UnauthorizedException.class)
    public void testGrant() {
        initAPI("0003");

            api.login("write", "write");
        //private
        realm.addPermissions(ClothoRealm.ANONYMOUS_USER, WRITE.actions, util.getPrivate().getId());
        //public
        realm.addPermissions(ClothoRealm.ANONYMOUS_USER, WRITE.actions, util.getPrivate().getId());        
    }
}
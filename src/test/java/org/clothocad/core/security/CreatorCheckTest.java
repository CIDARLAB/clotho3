/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import java.util.HashMap;
import java.util.Map;
import org.apache.shiro.SecurityUtils;
import org.apache.shiro.authz.*;

import org.apache.shiro.subject.Subject;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.model.Institution;
import org.junit.Ignore;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 * test case of user with 'write' permission test actions of read, edit, delete
 * and edit permission on a private object
 *
 * @author yu
 * @version 1.0
 */
public class CreatorCheckTest extends AbstractSecurityTest {

    /**
     * constructor
     */
    public CreatorCheckTest() {
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
    @Ignore @Test(expected = UnauthorizedException.class)
    public void testEditPermission() {
        initAPI("0003");

            Subject currentUser = SecurityUtils.getSubject();
            api.login("write", "write");
            /*
             *code here to edit permission 
             */

    }
}
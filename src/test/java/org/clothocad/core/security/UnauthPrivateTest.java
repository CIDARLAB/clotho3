/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import com.google.common.collect.Lists;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import javax.persistence.EntityNotFoundException;
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
 * test case of user with 'run' permission test actions of read, edit, delete
 * and edit permission on a private object
 *
 * @author yu
 * @version 1.0
 */
public class UnauthPrivateTest extends AbstractSecurityTest {

    /**
     * constructor
     */
    public UnauthPrivateTest() {
        super();
    }
    
    /**
     * test read action
     *
     */
    @Test(expected = EntityNotFoundException.class)
    public void testRead() {
        initAPI("0000");

        ObjBase priv = util.getPrivate();
        persistor.get(Institution.class, priv.getId());
    }

    @Test
    public void testFind() {
        initAPI("anonFind");
        Map<String, Object> query = new HashMap<>();
        query.put("schema", Institution.class.getName());
        assertEquals(1, Lists.newArrayList(persistor.find(query)).size());
    }
    
    @Test(expected = EntityNotFoundException.class)
    public void testRun() throws Exception {
        initAPI("anonRun");
        Map<String,Object> command = new HashMap<>();
        command.put("args", new ArrayList());
        command.put("id", util.getPrivateModule().getId());
        command.put("function", "function");
        api.run(command);
    }
    
    /**
     * test edit action
     *
     * @exception UnauthorizedException
     */
    @Test(expected = UnauthorizedException.class)
    public void testEdit() {
        initAPI("0001");
        ObjBase priv = util.getPrivate();
        try{

            Map<String,Object> edit = new HashMap<>();
            edit.put("id", priv.getId());
            edit.put("name", "Changed Name");
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
            /*
             *code here to edit permission 
             */

    }
}
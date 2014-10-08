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
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.model.Institution;
import org.junit.Test;
import static org.junit.Assert.*;
import static org.clothocad.core.security.ClothoPermission.*;

/**
 * test case of user with 'owner' permission test actions of read, edit, delete
 * and edit permission on a private object
 *
 * @author yu
 * @version 1.0
 * 
 * TODO: both ObjBase & JSON paths tested
 */
public class OwnPermissionTest extends AbstractSecurityTest {

    /**
     * constructor
     */
    public OwnPermissionTest() {
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
        api.login("owner", "owner");
        assertEquals(priv.getId(), persistor.get(Institution.class, priv.getId()).getId());
    }


    @Test
    public void testFind() {
        initAPI("ownerFind");
            api.login("owner", "owner");
        Map<String, Object> query = new HashMap<>();
        query.put("schema", Institution.class.getName());
        assertEquals(2, Lists.newArrayList(persistor.find(query)).size());
        }

    @Test
    public void testRun() throws Exception {
        initAPI("ownerRun");
        login("owner");
        Map<String,Object> command = new HashMap<>();
        command.put("args", new ArrayList());
        command.put("id", util.getPrivateModule().getId());
        command.put("function", "function");
        assertEquals("function ran!", api.run(command));
    }

    //TODO: test run script
    /**
     * test edit action
     */
    @Test
    public void testEdit() {
        initAPI("0001");
        ObjBase priv = util.getPrivate();
        try {

            Map<String, Object> edit = new HashMap<>();
            edit.put("id", priv.getId());
            edit.put("name", "Changed Name");
            api.login("owner", "owner");
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
     */
    @Test
    public void testDelete() {

        try {
                    initAPI("0002");
        ObjectId id = util.getPrivate().getId();
            api.login("owner", "owner");
            persistor.delete(id);
            try {
                persistor.getAsJSON(id);
                fail();
            } catch (EntityNotFoundException e){
                
            }

        } finally {
            persistor.save(util.getPrivate());
        }

    }

    /**
     * test edit permission
     */
    @Test
    public void testGrant() {
        initAPI("0003");

            api.login("owner", "owner");
        //private
        realm.addPermissions(ClothoRealm.ANONYMOUS_USER, WRITE.actions, util.getPrivate().getId());
        //public
        realm.addPermissions(ClothoRealm.ANONYMOUS_USER, WRITE.actions, util.getPrivate().getId());        
        }

    //TODO: check behavior when granting/removing public status
    //TODO: check behavior when granting redundant permissions
    //TODO: check for ownerless objects
}
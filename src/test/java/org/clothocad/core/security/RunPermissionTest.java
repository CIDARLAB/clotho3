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
public class RunPermissionTest extends AbstractSecurityTest {

    /**
     * constructor
     */
    public RunPermissionTest() {
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
        api.login("run", "run");
        persistor.get(Institution.class, priv.getId());
    }
    
    @Test
    public void testSneakyRead() throws Exception  {
        initAPI("sneakyRead");
        api.login("run", "run");        
        Map<String,Object> function = new HashMap<>();
        function.put("schema", "org.clothocad.core.datums.Function");
        function.put("code", " function (){\n" +
        "     return clotho.get(\""+ util.getPrivateModule().getId().toString() +"\");\n" +
        " }");
        function.put("language", "JAVASCRIPT");
        
        ObjectId functionId = api.create(function);
        Map<String,Object> arguments = new HashMap();
        arguments.put("id", functionId.toString());
        arguments.put("args", new String[0]);
        assertEquals(Void.TYPE, api.run(arguments));
    }

    @Test
    public void testFind() {
        initAPI("runFind");
        api.login("run", "run");
        Map<String, Object> query = new HashMap<>();
        query.put("schema", Institution.class.getName());
        assertEquals(1, Lists.newArrayList(persistor.find(query)).size());
    }
    
    @Test
    public void testRun() throws Exception {
        initAPI("runRun");
        login("run");
        Map<String,Object> command = new HashMap<>();
        command.put("args", new ArrayList());
        command.put("id", util.getPrivateModule().getId());
        command.put("function", "function");
        assertEquals("function ran!", api.run(command));
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
            api.login("run", "run");
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
            api.login("run", "run");
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
            api.login("run", "run");
            /*
             *code here to edit permission 
             */

    }
}
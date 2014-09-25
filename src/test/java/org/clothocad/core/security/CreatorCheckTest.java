/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import javax.persistence.EntityNotFoundException;
import org.apache.shiro.SecurityUtils;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.model.Institution;
import org.junit.Test;

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
    public void testCreatorIsOwner() {
        initAPI("0000");
        api.login("read", "read");
        ObjectId id = persistor.save(new Institution("New Institution", "", "", ""));
        SecurityUtils.getSubject().checkPermission("data:delete:"+id.toString());
        SecurityUtils.getSubject().checkPermission("data:grant:"+id.toString());
    }

    /**
     * test edit action
     *
     */
    @Test(expected = EntityNotFoundException.class)
    public void testCreatedIsPrivate() {
        initAPI("0000");
        api.login("read", "read");
        ObjectId id = persistor.save(new Institution("New Institution", "", "", ""));
        api.logout();
        api.login("owner", "owner");
        persistor.get(Institution.class, id);
    }

}
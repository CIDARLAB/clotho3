/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.communication;

import com.google.common.collect.Sets;
import java.util.HashSet;
import java.util.Set;
import org.apache.shiro.authz.AuthorizationException;
import org.apache.shiro.subject.PrincipalCollection;
import org.apache.shiro.subject.SimplePrincipalCollection;
import static org.clothocad.core.communication.AbstractServerAPITest.api;
import static org.clothocad.core.communication.AbstractServerAPITest.realm;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.security.ClothoAction;
import static org.junit.Assert.*;
import org.junit.Test;
import static org.clothocad.core.security.ClothoAction.*;
import org.clothocad.core.security.ClothoPermission;
import org.clothocad.core.security.ClothoRealm;
import static org.clothocad.core.security.PermissionHolder.constructPermissionString;
import org.clothocad.model.Institution;

/**
 *
 * @author spaige
 */
public class GrantTest extends AbstractServerAPITest {

    @Override
    public void setUp() {
        super.setUp();
    }

    @Test
    public void testAddPublic() {
        ObjectId id = makeInstitution();
        api.grant(id, null, Sets.newHashSet("public"), new HashSet<String>());
        assertHasOnlyPermissions(id, "none", view, run);
    }

    @Test
    public void testRemovePublic() {
        ObjectId id = makeInstitution();
        realm.setPublic(id);
        api.grant(id, null, new HashSet<String>(), Sets.newHashSet("public"));
        assertHasOnlyPermissions(id, "none");
    }

    @Test
    public void testAddRead() {
        ObjectId id = makeInstitution();
        api.grant(id, "none", Sets.newHashSet("read"), new HashSet<String>());
        assertHasOnlyPermissions(id, "none", view, run);
    }

    @Test
    public void testRemoveRead() {
        ObjectId id = makeInstitution();
        realm.addPermissions("none", ClothoPermission.OWN.actions, id);
        api.grant(id, "none", new HashSet<String>(), Sets.newHashSet("read"));
        assertHasOnlyPermissions(id, "none", run);
    }

    @Test
    public void testAddEdit() {
         ObjectId id = makeInstitution();
        api.grant(id, "none", Sets.newHashSet("write"), new HashSet<String>());
        assertHasOnlyPermissions(id, "none", run, view, edit);
    }

    @Test
    public void testRemoveEdit() {
        ObjectId id = makeInstitution();
        realm.addPermissions("none", ClothoPermission.OWN.actions, id);
        api.grant(id, "none", new HashSet<String>(), Sets.newHashSet("write"));
        assertHasOnlyPermissions(id, "none", view, run);
    }

    @Test
    public void testAddRun() {
        ObjectId id = makeInstitution();
        api.grant(id, "none", Sets.newHashSet("run"), new HashSet<String>());
        assertHasOnlyPermissions(id, "none", run);
    }

    @Test
    public void testRemoveRun() {
        ObjectId id = makeInstitution();
        realm.addPermissions("none", ClothoPermission.OWN.actions, id);
        api.grant(id, "none", new HashSet<String>(), Sets.newHashSet("run"));
        assertHasOnlyPermissions(id, "none");
    }

    @Test
    public void testAddOwn() {
        ObjectId id = makeInstitution();
        api.grant(id, "none", Sets.newHashSet("own"), new HashSet<String>());
        assertHasOnlyPermissions(id, "none", view, run, edit, delete, grant);
    }

    @Test
    public void testRemoveOwn() {
        ObjectId id = makeInstitution();
        realm.addPermissions("none", ClothoPermission.OWN.actions, id);
        api.grant(id, "none", new HashSet<String>(), Sets.newHashSet("own"));
        assertHasOnlyPermissions(id, "none", view, run, edit);
    }

    private ObjectId makeInstitution() {
        return persistor.save(new Institution("Test Institution", "", "", ""));
    }

    private void assertHasOnlyPermissions(ObjectId id, String principal, ClothoAction... clothoActions){
        assertHasOnlyPermissions(id, principal, realm, clothoActions);
    }
    
    public static void assertHasOnlyPermissions(ObjectId id, String principal, ClothoRealm realm, ClothoAction... clothoActions) {
        Set<ClothoAction> shouldHavePermissions = Sets.newHashSet(clothoActions);
        Set<ClothoAction> shouldLackPermissions = Sets.complementOf(shouldHavePermissions, ClothoAction.class);
        Set<ClothoAction> hasPermissions = new HashSet<>();
        Set<ClothoAction> lacksPermissions = new HashSet<>();
        PrincipalCollection principals = new SimplePrincipalCollection(principal, "clotho");
        for (ClothoAction action : ClothoAction.values()) {
            try {
                realm.checkPermission(principals, constructPermissionString(action, id));
                hasPermissions.add(action);
            } catch (AuthorizationException e) {
                lacksPermissions.add(action);
            }
        }
        assertEquals(shouldHavePermissions, hasPermissions);
        assertEquals(shouldLackPermissions, lacksPermissions);
    }
}

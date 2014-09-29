/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.communication;

import com.google.common.collect.Sets;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.security.ClothoAction;
import static org.clothocad.core.security.ClothoAction.delete;
import static org.clothocad.core.security.ClothoAction.edit;
import static org.clothocad.core.security.ClothoAction.grant;
import static org.clothocad.core.security.ClothoAction.run;
import static org.clothocad.core.security.ClothoAction.view;
import org.clothocad.core.security.ClothoPermission;
import org.clothocad.core.security.ClothoRealm;
import org.clothocad.model.Institution;
import org.junit.Test;

/**
 *
 * @author spaige
 */
public class RouterGrantTest extends AbstractRouterTest {

    protected ClothoRealm realm;
    protected Persistor persistor;
    
    public RouterGrantTest() {
        super();
        realm = injector.getInstance(ClothoRealm.class);
        persistor = injector.getInstance(Persistor.class);
    }
    
    @Test
    public void testAddPublic() throws IOException{
        ObjectId id = makeInstitution();
        sendGrantMessage(id, null, Sets.newHashSet("public"), new HashSet<String>());
        assertHasOnlyPermissions(id, "none", view, run);
    }

    @Test
    public void testRemovePublic() throws IOException{
        ObjectId id = makeInstitution();
        realm.setPublic(id);
        sendGrantMessage(id, null, new HashSet<String>(), Sets.newHashSet("public"));
        assertHasOnlyPermissions(id, "none");
    }

    @Test
    public void testAddRead() throws IOException{
        ObjectId id = makeInstitution();
        sendGrantMessage(id, "none", Sets.newHashSet("read"), new HashSet<String>());
        assertHasOnlyPermissions(id, "none", view, run);
    }

    @Test
    public void testRemoveRead() throws IOException{
        ObjectId id = makeInstitution();
        realm.addPermissions("none", ClothoPermission.OWN.actions, id);
        sendGrantMessage(id, "none", new HashSet<String>(), Sets.newHashSet("read"));
        assertHasOnlyPermissions(id, "none", run);
    }

    @Test
    public void testAddEdit() throws IOException{
         ObjectId id = makeInstitution();
        sendGrantMessage(id, "none", Sets.newHashSet("write"), new HashSet<String>());
        assertHasOnlyPermissions(id, "none", run, view, edit);
    }

    @Test
    public void testRemoveEdit() throws IOException{
        ObjectId id = makeInstitution();
        realm.addPermissions("none", ClothoPermission.OWN.actions, id);
        sendGrantMessage(id, "none", new HashSet<String>(), Sets.newHashSet("write"));
        assertHasOnlyPermissions(id, "none", view, run);
    }

    @Test
    public void testAddRun() throws IOException{
        ObjectId id = makeInstitution();
        sendGrantMessage(id, "none", Sets.newHashSet("run"), new HashSet<String>());
        assertHasOnlyPermissions(id, "none", run);
    }

    @Test
    public void testRemoveRun() throws IOException{
        ObjectId id = makeInstitution();
        realm.addPermissions("none", ClothoPermission.OWN.actions, id);
        sendGrantMessage(id, "none", new HashSet<String>(), Sets.newHashSet("run"));
        assertHasOnlyPermissions(id, "none");
    }

    @Test
    public void testAddOwn() throws IOException{
        ObjectId id = makeInstitution();
        sendGrantMessage(id, "none", Sets.newHashSet("own"), new HashSet<String>());
        assertHasOnlyPermissions(id, "none", view, run, edit, delete, grant);
    }

    @Test
    public void testRemoveOwn() throws IOException{
        ObjectId id = makeInstitution();
        realm.addPermissions("none", ClothoPermission.OWN.actions, id);
        sendGrantMessage(id, "none", new HashSet<String>(), Sets.newHashSet("own"));
        assertHasOnlyPermissions(id, "none", view, run, edit);
    }

    
    private void assertHasOnlyPermissions(ObjectId id, String principal, ClothoAction... clothoActions){
        GrantTest.assertHasOnlyPermissions(id, principal, realm, clothoActions);
    }
    
    private ObjectId makeInstitution() {
        return persistor.save(new Institution("Test Institution", "", "", ""));
    }

    private void sendGrantMessage(ObjectId id, String user, Set<String> add, Set<String> remove) throws IOException {
        TestConnection connection = new TestConnection("");
        Map<String, Object> data = new HashMap<>();
        data.put("id", id.toString());
        data.put("add", add);
        data.put("remove", remove);
        data.put("user", user);
        
        sendMessage(new Message(Channel.grant, data, ""),connection);
    }
    
}

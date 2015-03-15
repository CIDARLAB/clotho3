/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import org.apache.shiro.SecurityUtils;
import org.apache.shiro.subject.SimplePrincipalCollection;
import org.clothocad.core.datums.Module;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.security.ClothoRealm;
import org.clothocad.core.security.Visibility;
import org.clothocad.model.Institution;
import static org.clothocad.core.security.ClothoPermission.*;
import org.clothocad.model.Person;

/**
 *
 * @author spaige
 *
 * Provides a simple universe of objects and users for permission testing
 */
public class SecurityTestUtils {

    public SecurityTestUtils() {
    }
    /*
     Holds credentials for 4 users in <username, password> format. 
     usernames are 'none', 'read', 'write', and 'owner', and describe the permission levels
     each user has on both the private and public object in the universe.
     */
    public Map<String, String> credentials;
    public Map<String, ObjBase> objects;

    public ObjBase getPublic() {
        return objects.get("public");
    }

    public ObjBase getPrivate() {
        return objects.get("private");
    }

    public Module getPublicModule() {
        return (Module) objects.get("publicModule");
    }

    public Module getPrivateModule() {
        return (Module) objects.get("privateModule");
    }

    private void doCreateTestRealmData(Persistor p, ClothoRealm realm) {
        //make objects
        Institution privateInstitution = new Institution("Private Institution", "", "", "");
        privateInstitution.setVisibility(Visibility.PRIVATE);
        p.save(privateInstitution);
        Institution publicInstitution = new Institution("Public Institution", "", "", "");
        p.save(publicInstitution);
        //XXX: need dummy function
        Module privateModule = new SecurityTester();
        //privateModule.setId(new ObjectId(privateModule.getId().toString() + "$Private"));
        p.save(privateModule);
        Module publicModule = new SecurityTester();
        //publicModule.setId(new ObjectId(publicModule.getId().toString() + "$Public"));
        p.save(publicModule);

        objects = new HashMap<>();
        objects.put(("private"), privateInstitution);
        objects.put(("public"), publicInstitution);
        objects.put(("privateModule"), privateModule);
        objects.put(("publicModule"), publicModule);

        Set<ObjectId> ids = new HashSet<>();
        ids.add(privateInstitution.getId());
        ids.add(publicInstitution.getId());
        ids.add(publicModule.getId());
        ids.add(privateModule.getId());

        //make users
        createAccount(p, realm, "none");
        createAccount(p, realm, "read");
        realm.addPermissions("read", READ.actions, ids);
        createAccount(p, realm, "write");
        realm.addPermissions("write", WRITE.actions, ids);
        createAccount(p, realm, "run");
        realm.addPermissions("run", RUN.actions, ids);
        createAccount(p, realm, "owner");
        realm.addPermissions("owner", OWN.actions, ids);

        //make public objects public
        realm.setPublic(publicInstitution.getId());
        realm.setPublic(publicModule.getId());

        credentials = new HashMap<>();
        credentials.put("none", "none");
        credentials.put("read", "read");
        credentials.put("write", "write");
        credentials.put("run", "run");
        credentials.put("owner", "owner");

    }

    public static class CreateTestRealmData implements Callable<SecurityTestUtils> {

        Persistor p;
        ClothoRealm realm;

        public CreateTestRealmData(Persistor p, ClothoRealm realm) {
            this.p = p;
            this.realm = realm;
        }

        @Override
        public SecurityTestUtils call() throws Exception {
            SecurityTestUtils util = new SecurityTestUtils();
            util.doCreateTestRealmData(p, realm);
            return util;
        }
    }

    public static void createAccount(Persistor p, ClothoRealm realm, String account) {

        realm.addAccount(account, account);

        Person person = new Person(account + "-permission User");
        person.setPrimaryEmail(account);
        person.setEmailAddress(account);
        person.setId(new ObjectId());
        ObjectId personId = p.save(person);

        SecurityUtils.getSubject().runAs(new SimplePrincipalCollection(account, "clotho"));
        realm.addPrincipal(account, personId, "persons");
        SecurityUtils.getSubject().releaseRunAs();
        
        realm.addPermissions(account, OWN.actions, personId);
    }
}

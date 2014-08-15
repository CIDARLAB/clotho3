/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import com.google.common.collect.ImmutableSet;
import java.util.List;
import java.util.Set;
import javax.inject.Inject;
import lombok.extern.slf4j.Slf4j;
import org.apache.shiro.authc.AccountException;
import org.apache.shiro.authc.AuthenticationException;
import org.apache.shiro.authc.AuthenticationInfo;
import org.apache.shiro.authc.AuthenticationToken;
import org.apache.shiro.authc.UsernamePasswordToken;
import org.apache.shiro.authc.credential.HashedCredentialsMatcher;
import org.apache.shiro.authz.AuthorizationInfo;
import org.apache.shiro.authz.permission.RolePermissionResolver;
import org.apache.shiro.crypto.hash.Sha256Hash;
import org.apache.shiro.realm.AuthorizingRealm;
import org.apache.shiro.subject.PrincipalCollection;
import org.apache.shiro.util.CollectionUtils;
import org.clothocad.core.datums.ObjectId;
import static org.clothocad.core.security.ServerSubject.SERVER_USER;

/**
 *
 * Note: probably worth implementing a more intelligent permission search for
 * data object permissions
 * 
 * TODO: locking
 *
        //TODO: Networking: realm name is clotho instance name
        //                  add new user to local instance group 
 * @author spaige
 */
@Slf4j
public class ClothoRealm extends AuthorizingRealm {

    public static final Set<String> READ = ImmutableSet.of("view", "run");
    public static final Set<String> WRITE = ImmutableSet.of("view", "edit", "run");
    public static final Set<String> RUN = ImmutableSet.of("run");
    public static final Set<String> OWN = ImmutableSet.of("view", "run", "edit", "delete", "grant");
    
    public static final String ALL = "_all";
    
    public static final String ANONYMOUS_USER = "_anonymous";
    
    private CredentialStore store;
    
    public static final AuthenticationToken getAnonymousUserToken(){
        return new UsernamePasswordToken(ANONYMOUS_USER, ANONYMOUS_USER);
    }

    @Inject
    public ClothoRealm(CredentialStore store, RolePermissionResolver roleResolver) {
        super();

        //XXX: up number of iterations
        HashedCredentialsMatcher matcher = new HashedCredentialsMatcher(Sha256Hash.ALGORITHM_NAME);
        matcher.setStoredCredentialsHexEncoded(false);

        this.store = store;
        setAuthenticationTokenClass(UsernamePasswordToken.class);
        setCredentialsMatcher(matcher);

        setRolePermissionResolver(roleResolver);

        setUpRealm();
    }
    
    protected void setUpRealm(){
        //set up default groups
        if (store.getGroup(ALL) == null){
            addGroup(ALL);
        }       
        
        //set up anonymous user
        if (store.getAccount(ANONYMOUS_USER) == null){
            addAccount(ANONYMOUS_USER, ANONYMOUS_USER);
        }
        
        if (store.getAccount(SERVER_USER) == null){
            store.saveAccount(new DummyAccount(SERVER_USER));
        }
    }

    @Override
    protected AuthorizationInfo doGetAuthorizationInfo(PrincipalCollection pc) {
        log.debug("getting authz info for {}", pc);

        ClothoAccount account = store.getAccount(pc.getPrimaryPrincipal().toString());
        return account.getAuthzInfo();
    }

    @Override
    protected AuthenticationInfo doGetAuthenticationInfo(AuthenticationToken at) throws AuthenticationException {
        log.debug("getting authc info for {}", at);

        ClothoAccount account = store.getAccount(((UsernamePasswordToken) at).getUsername());
        if (!account.isAuthenticatable()) throw new AccountException("Cannot authenticate as " + at.getPrincipal().toString());
        return account;
    }

    public void addAccount(String username, String password) {
        if (store.getAccount(username) != null){
            throw new javax.persistence.EntityExistsException();
        }
        store.saveAccount(new ClothoAccount(username, password));
    }

    
    public void addGroup(String groupName){
        store.saveGroup(new AuthGroup(groupName));
    }
    

    public void addPermission(String username, String permission) {
        ClothoAccount account = store.getAccount(username);
        //parts list: data, <object id>, <permission>
        account.getAuthzInfo().addPermission(permission);
        store.saveAccount(account);
    }

    public void addPermissionToGroup(String groupName, String permission) {
        AuthGroup group = store.getGroup(groupName);
        group.addPermission(permission);
        store.saveGroup(group);
    }

    public void removePermissionFromGroup(String groupName, String permission) {
        AuthGroup group = store.getGroup(groupName);
        group.removePermission(permission);
        store.saveGroup(group);
    }

    public void removePermission(String username, String permission) {
        ClothoAccount account = store.getAccount(username);
        //parts list: data, <object id>, <permission>
        account.getAuthzInfo().removePermission(permission);
        store.saveAccount(account);
    }

    public void addPermissions(String username, Set<String> permissions, Set<ObjectId> ids) {
        for (String permission : permissions) {
            for (ObjectId id : ids) {
                addPermission(username, "data:" + permission + ":" + id.toString());
            }
        }
    }

    public void addPermissions(String username, Set<String> permissions, ObjectId id) {
        for (String permission : permissions) {
            addPermission(username, "data:" + permission + ":" + id.toString());
        }
    }

    public void addPermissionsToGroup(String groupName, Set<String> permissions, ObjectId id){
        for (String permission : permissions){
            addPermissionToGroup(groupName, "data:" + permission +":" +id.toString());
        }
    }
    
    public void removePermissionsFromGroup(String groupName, Set<String> permissions, ObjectId id){
        for (String permission : permissions){
            removePermissionFromGroup(groupName, "data:" + permission + ":" + id.toString());
        }
    }
    
    public void deleteAll() {
        store.deleteAllCredentials();
        setUpRealm();
    }

    public void removePublic(ObjectId id) {
        removePermissionsFromGroup(ALL,READ,id);
    }

    public void setPublic(ObjectId id) {
        addPermissionsToGroup(ALL,READ,id);
    }
    
    public String getRealmName() {
        return "clotho";
    }
}

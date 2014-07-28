/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import com.google.common.collect.ImmutableSet;
import java.util.HashSet;
import java.util.Set;
import javax.inject.Inject;
import lombok.extern.slf4j.Slf4j;
import org.apache.shiro.authc.Account;
import org.apache.shiro.authc.AuthenticationException;
import org.apache.shiro.authc.AuthenticationInfo;
import org.apache.shiro.authc.AuthenticationToken;
import org.apache.shiro.authc.UsernamePasswordToken;
import org.apache.shiro.authc.credential.HashedCredentialsMatcher;
import org.apache.shiro.authz.AuthorizationInfo;
import org.apache.shiro.authz.SimpleAuthorizationInfo;
import org.apache.shiro.authz.permission.RolePermissionResolver;
import org.apache.shiro.crypto.SecureRandomNumberGenerator;
import org.apache.shiro.crypto.hash.Sha256Hash;
import org.apache.shiro.crypto.hash.SimpleHash;
import org.apache.shiro.realm.AuthorizingRealm;
import org.apache.shiro.subject.PrincipalCollection;
import org.apache.shiro.util.ByteSource;
import org.clothocad.core.datums.ObjectId;

/**
 *
 * Note: probably worth implementing a more intelligent permission search for 
 * data object permissions
 * 
 * @author spaige
 */
@Slf4j
public class ClothoRealm extends AuthorizingRealm {
        
    public static final Set<String> READ = ImmutableSet.of("view", "run");
    public static final Set<String> WRITE = ImmutableSet.of("view", "edit", "run");
    public static final Set<String> RUN = ImmutableSet.of("run");
    public static final Set<String> OWN = ImmutableSet.of("view", "run", "edit", "delete", "grant");
    
    private CredentialStore store;
    
    @Inject
    public ClothoRealm (CredentialStore store, RolePermissionResolver roleResolver){
        super();
        
        //XXX: up number of iterations
        HashedCredentialsMatcher matcher = new HashedCredentialsMatcher(Sha256Hash.ALGORITHM_NAME);
        matcher.setStoredCredentialsHexEncoded(false);
        
        this.store = store;
        setAuthenticationTokenClass(UsernamePasswordToken.class);
        setCredentialsMatcher(matcher);
        
        setRolePermissionResolver(roleResolver);
        
        //set up role info
    }

    @Override
    protected AuthorizationInfo doGetAuthorizationInfo(PrincipalCollection pc) {
        log.debug("getting authz info for {}", pc);
        Account account = store.getAccount(pc.getPrimaryPrincipal().toString());
        SimpleAuthorizationInfo authInfo =  new SimpleAuthorizationInfo();
        //authInfo.setObjectPermissions(new HashSet<>(account.getObjectPermissions()));
        authInfo.setRoles(new HashSet<>(account.getRoles()));
        authInfo.addRole("any");
        authInfo.setStringPermissions(new HashSet<>(account.getStringPermissions()));
        return authInfo;
    }

    @Override
    protected AuthenticationInfo doGetAuthenticationInfo(AuthenticationToken at) throws AuthenticationException {
        log.debug("getting authc info for {}", at);
        
        return store.getAccount(((UsernamePasswordToken) at).getUsername());
    }
    
    //XXX: check for preexisting account
    public void addAccount(String username, String password) {
        
        ByteSource salt = new SecureRandomNumberGenerator().nextBytes();
        SimpleHash hashedPw = new SimpleHash(Sha256Hash.ALGORITHM_NAME, password, salt);
        
        store.saveAccount(username, hashedPw, salt);
    }
    
    public void addPermission(String username, String permission){
        store.addPermission(username, permission);
    }

    public void removePermission(String username, String permission){
        store.removePermission(username, permission);
    }
    
    
    public void addPermissions(String username, Set<String> permissions, Set<ObjectId> ids){
        for (String permission : permissions){
            for (ObjectId id : ids){
                addPermission(username, "data:" + permission + ":" + id.toString());
            }
        }
    }

    public void addPermissions(String username, Set<String> permissions, ObjectId id){
        for (String permission : permissions){
            addPermission(username, "data:" + permission + ":" + id.toString());
        }
    }
    
    public void deleteAll(){
        store.deleteAllCredentials();
    }
}

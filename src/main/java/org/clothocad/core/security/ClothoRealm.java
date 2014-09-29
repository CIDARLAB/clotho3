/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import javax.inject.Inject;
import lombok.extern.slf4j.Slf4j;
import org.apache.shiro.authc.AuthenticationException;
import org.apache.shiro.authc.AuthenticationInfo;
import org.apache.shiro.authc.AuthenticationToken;
import org.apache.shiro.authc.SimpleAccount;
import org.apache.shiro.authc.UsernamePasswordToken;
import org.apache.shiro.authc.credential.CredentialsMatcher;
import org.apache.shiro.authc.credential.HashedCredentialsMatcher;
import org.apache.shiro.authz.AuthorizationInfo;
import org.apache.shiro.crypto.SecureRandomNumberGenerator;
import org.apache.shiro.crypto.hash.Sha256Hash;
import org.apache.shiro.crypto.hash.SimpleHash;
import org.apache.shiro.realm.AuthorizingRealm;
import org.apache.shiro.subject.PrincipalCollection;
import org.apache.shiro.util.ByteSource;
import org.clothocad.core.persistence.Persistor;

/**
 *
 * @author spaige
 */
@Slf4j
public class ClothoRealm extends AuthorizingRealm {
    
    private CredentialStore store;
    
    @Inject
    public ClothoRealm (CredentialStore store){
        super();
        
        //XXX: up number of iterations
        HashedCredentialsMatcher matcher = new HashedCredentialsMatcher(Sha256Hash.ALGORITHM_NAME);
        matcher.setStoredCredentialsHexEncoded(false);
        
        this.store = store;
        setAuthenticationTokenClass(UsernamePasswordToken.class);
        setCredentialsMatcher(matcher);
    }

    @Override
    protected AuthorizationInfo doGetAuthorizationInfo(PrincipalCollection pc) {
        log.debug("getting authz info for {}", pc);
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    protected AuthenticationInfo doGetAuthenticationInfo(AuthenticationToken at) throws AuthenticationException {
        log.debug("getting authc info for {}", at);
        
        return store.getAccount(((UsernamePasswordToken) at).getUsername());
    }
    
    public void addAccount(String username, String password) {
        
        ByteSource salt = new SecureRandomNumberGenerator().nextBytes();
        SimpleHash hashedPw = new SimpleHash(Sha256Hash.ALGORITHM_NAME, password, salt);
        
        store.saveAccount(username, hashedPw, salt);
        
    }

    public void deleteAll(){
        store.deleteAllCredentials();
    }
    public void updatePassword(String username, String password)
    {
        ByteSource salt = new SecureRandomNumberGenerator().nextBytes();
        SimpleHash hashedPw = new SimpleHash(Sha256Hash.ALGORITHM_NAME, password, salt);
        System.out.println("Password "+store.getAccount(username).getCredentials().toString());
         
        store.updatePassword(username, hashedPw, salt);
       
    }
}

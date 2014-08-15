/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;


import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonTypeInfo;
import lombok.Delegate;
import lombok.Getter;
import org.apache.shiro.authc.Account;
import org.apache.shiro.authc.MergableAuthenticationInfo;
import org.apache.shiro.authc.SaltedAuthenticationInfo;
import org.apache.shiro.authc.SimpleAuthenticationInfo;
import org.apache.shiro.crypto.SecureRandomNumberGenerator;
import org.apache.shiro.crypto.hash.Sha256Hash;
import org.apache.shiro.crypto.hash.SimpleHash;
import org.apache.shiro.subject.SimplePrincipalCollection;
import org.apache.shiro.util.ByteSource;
import org.jongo.marshall.jackson.oid.Id;

/**
 *
 * @author spaige
 */
@JsonTypeInfo(use=JsonTypeInfo.Id.CLASS, include=JsonTypeInfo.As.PROPERTY, property="@class")
public class ClothoAccount implements Account, MergableAuthenticationInfo, SaltedAuthenticationInfo {

    
    @JsonCreator
    public ClothoAccount(@JsonProperty("_id") String username,
    @JsonProperty("authcInfo") SimpleAuthenticationInfo authcInfo, 
    @JsonProperty("authzInfo") ClothoAuthorizationInfo authzInfo) {
        this.id = username;
        this.authcInfo = authcInfo;
        this.authzInfo = authzInfo;
    }
    
    public boolean isAuthenticatable(){
        return true;
    }
    
    
    
    public ClothoAccount(String username, String password){
        this(username);
            
        ByteSource salt = new SecureRandomNumberGenerator().nextBytes();
        SimpleHash hashedPw = new SimpleHash(Sha256Hash.ALGORITHM_NAME, password, salt);
        authcInfo = new SimpleAuthenticationInfo(new SimplePrincipalCollection(username,"clotho"), hashedPw, salt);

    }
    
    public ClothoAccount(String username){
        id = username;
        authzInfo = new ClothoAuthorizationInfo();
        authzInfo.addGroup(ClothoRealm.ALL);      
    }    
    
    @Id
    private String id;
    
    @Delegate
    private SimpleAuthenticationInfo authcInfo;
    
    @Delegate
    @Getter
    private ClothoAuthorizationInfo authzInfo;

}

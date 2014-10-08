/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import org.apache.shiro.crypto.hash.SimpleHash;
import org.apache.shiro.util.ByteSource;

/**
 *
 * @author spaige
 */
public interface CredentialStore {
    
    public ClothoAccount getAccount(String username);
      
    public void saveAccount(ClothoAccount account);

    public AuthGroup getGroup(String groupName);        
    public void deleteAllCredentials();

    public void saveGroup(AuthGroup authGroup);
}

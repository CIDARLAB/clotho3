/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import org.apache.shiro.authc.SimpleAccount;
import org.apache.shiro.crypto.hash.SimpleHash;
import org.apache.shiro.util.ByteSource;

/**
 *
 * @author spaige
 */
public interface CredentialStore {
    SimpleAccount getAccount(String username);
    
    public void saveAccount(String username, SimpleHash hashedPw, ByteSource salt);
    
    public void addPermission(String username, String permission);
    
    public void removePermission(String username, String permission);
    
    public void deleteAllCredentials();
}

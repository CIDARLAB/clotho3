/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import org.apache.shiro.authc.SimpleAccount;

/**
 *
 * @author spaige
 */
public interface CredentialStore {
    SimpleAccount getAccount(String username);
    
    void saveAccount(SimpleAccount account); 
}

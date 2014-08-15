/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

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

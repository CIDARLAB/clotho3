/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.clothocad.core.util;

import java.util.HashMap;
import java.util.Map;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.security.Visibility;
import org.clothocad.model.Institution;

/**
 *
 * @author spaige
 * 
 * Provides a simple universe of objects and users for permission testing
 */
public class SecurityTestUtils {
    
    public SecurityTestUtils(Persistor p){
        //make objects
        Institution privateInstitution = new Institution("Private Institution", "", "", "");
        privateInstitution.setVisibility(Visibility.PRIVATE);
        p.save(privateInstitution);
        Institution publicInstitution = new Institution("Public Institution", "", "", "");
        privateInstitution.setVisibility(Visibility.PUBLIC);
        p.save(publicInstitution);
        objects = new HashMap<>();
        objects.put(("private"), privateInstitution);
        objects.put(("public"), publicInstitution);
        
        //make users - just credentials for now
        credentials = new HashMap<>();
        credentials.put("none", "none");
        credentials.put("read", "read");
        credentials.put("write", "write");
        credentials.put("owner", "owner");
    }
    
    /*
    Holds credentials for 4 users in <username, password> format. 
    usernames are 'none', 'read', 'write', and 'owner', and describe the permission levels
    each user has on both the private and public object in the universe.
    */
    public Map<String,String>  credentials;
    public Map<String, ObjBase> objects;
}

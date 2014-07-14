/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.clothocad.core.util;

import com.google.common.collect.ImmutableSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.security.ClothoRealm;
import org.clothocad.core.security.Visibility;
import org.clothocad.model.Institution;

/**
 *
 * @author spaige
 * 
 * Provides a simple universe of objects and users for permission testing
 */
public class SecurityTestUtils {
    
    public SecurityTestUtils(Persistor p, ClothoRealm realm){
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
        
        Set<ObjectId> ids = new HashSet<>();
        ids.add(privateInstitution.getId());
        ids.add(publicInstitution.getId());
        
        //make users - just credentials for now
        realm.addAccount("none", "none");
        realm.addAccount("read", "read");
        addPermission("read", READ, ids, realm);
        realm.addAccount("write", "write");
        addPermission("write", WRITE, ids, realm);
        realm.addAccount("run", "run");
        addPermission("run", RUN, ids, realm);
        realm.addAccount("owner", "owner");
        addPermission("owner", OWN, ids, realm);
        
        credentials = new HashMap<>();
        credentials.put("none", "none");
        credentials.put("read", "read");
        credentials.put("write", "write");
        credentials.put("run", "run");
        credentials.put("owner", "owner");
    }
    
    private static void addPermission(String username, Set<String> permissions, Set<ObjectId> ids, ClothoRealm realm){
        for (String permission : permissions){
            for (ObjectId id : ids){
                realm.addPermission(username, "data:" + permission + ":" + id.toString());
            }
        }
    }
    
    /*
    Holds credentials for 4 users in <username, password> format. 
    usernames are 'none', 'read', 'write', and 'owner', and describe the permission levels
    each user has on both the private and public object in the universe.
    */
    public Map<String,String>  credentials;
    public Map<String, ObjBase> objects;
    
    public static final Set<String> READ = ImmutableSet.of("view", "run");
    public static final Set<String> WRITE = ImmutableSet.of("view", "edit", "run");
    public static final Set<String> RUN = ImmutableSet.of("run");
    public static final Set<String> OWN = ImmutableSet.of("view", "run", "edit", "delete", "grant");
}

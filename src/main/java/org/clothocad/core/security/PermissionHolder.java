/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import java.util.List;
import java.util.Map;
import lombok.Getter;
import org.apache.shiro.util.CollectionUtils;
import org.clothocad.core.datums.ObjectId;

/**
 *
 * @author spaige
 */
public abstract class PermissionHolder {
    
    @Getter
    protected Map<ObjectId, PermissionsOnObject> permissions;
    
    public static ObjectId getObjectId(String permission) {
        List<String> parts = CollectionUtils.asList(permission.split(":"));
        return new ObjectId(parts.get(2));
    }

    
    public void addPermission(String permission){
        ObjectId id = getObjectId(permission);
        if (!permissions.containsKey(id)){
            permissions.put(id, new PermissionsOnObject());
        }
        permissions.get(id).permissions.add(permission);
    }
    
    public void removePermission(String permission){
        ObjectId id = getObjectId(permission);
        if (permissions.containsKey(id)){
            permissions.get(id).permissions.remove(permission);
        }
    }    
}

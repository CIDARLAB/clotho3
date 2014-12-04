/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import org.apache.shiro.authz.AuthorizationException;
import org.apache.shiro.authz.Permission;
import org.apache.shiro.realm.Realm;
import org.apache.shiro.session.mgt.DefaultSessionManager;
import org.apache.shiro.subject.PrincipalCollection;
import org.apache.shiro.web.mgt.DefaultWebSecurityManager;
import org.clothocad.core.util.ClothoAuthoringEnvironment.AuthoringSubjectFactory;

/**
 *
 * @author spaige
 */
public class PermissiveSecurityManager extends DefaultWebSecurityManager {

    public PermissiveSecurityManager() {
        super();
        setSubjectFactory(new AuthoringSubjectFactory());
        setSessionManager(new DefaultSessionManager());
    }

    @SuppressWarnings({"UnusedDeclaration"})
    public PermissiveSecurityManager(Realm singleRealm) {
        this();
        setRealm(singleRealm);
    }

    @SuppressWarnings({"UnusedDeclaration"})
    public PermissiveSecurityManager(Collection<Realm> realms) {
        this();
        setRealms(realms);
    }

    @Override
    public void checkPermission(PrincipalCollection principals, Permission permission) throws AuthorizationException {
    }

    @Override
    public void checkPermission(PrincipalCollection principals, String permission) throws AuthorizationException {
    }

    @Override
    public void checkPermissions(PrincipalCollection principals, Collection<Permission> permissions) throws AuthorizationException {
    }

    @Override
    public void checkPermissions(PrincipalCollection principals, String... permissions) throws AuthorizationException {
    }

    @Override
    public boolean[] isPermitted(PrincipalCollection principals, List<Permission> permissions) {
        boolean[] trues = new boolean[permissions.size()];
        Arrays.fill(trues,true);
        return trues;
    }

    @Override
    public boolean isPermitted(PrincipalCollection principals, Permission permission) {
        return true;
    }

    @Override
    public boolean isPermitted(PrincipalCollection principals, String permissionString) {
        return true;
    }

    @Override
    public boolean[] isPermitted(PrincipalCollection principals, String... permissions) {
        boolean[] trues = new boolean[permissions.length];
        Arrays.fill(trues,true);
        return trues;       
    }

    @Override
    public boolean isPermittedAll(PrincipalCollection principals, Collection<Permission> permissions) {
        return true;
    }

    @Override
    public boolean isPermittedAll(PrincipalCollection principals, String... permissions) {
        return true;
    }

}
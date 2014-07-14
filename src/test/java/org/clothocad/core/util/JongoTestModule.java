/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import com.google.inject.AbstractModule;
import com.google.inject.Singleton;
import org.apache.shiro.authz.permission.RolePermissionResolver;
import org.clothocad.core.persistence.ClothoConnection;
import org.clothocad.core.security.CredentialStore;

/**
 *
 * @author spaige
 */
public class JongoTestModule extends AbstractModule {

    @Override
    protected void configure() {
        bind(CredentialStore.class).to(TestEnvConnection.class);
        bind(TestEnvConnection.class).in(Singleton.class);
        bind(ClothoConnection.class).to(TestEnvConnection.class);
        bind(RolePermissionResolver.class).to(TestEnvConnection.class);
    }

}

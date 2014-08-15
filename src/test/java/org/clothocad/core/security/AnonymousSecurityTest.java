/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import org.apache.shiro.SecurityUtils;
import org.apache.shiro.subject.SimplePrincipalCollection;
import org.apache.shiro.subject.Subject;
import static org.clothocad.core.security.AbstractShiroTest.getSecurityManager;

/**
 *
 * @author spaige
 */
public class AnonymousSecurityTest extends AbstractSecurityTest {
    public AnonymousSecurityTest() {
        super();  
    }
    
    @Override
    public void beforeTest(){
        Subject subjectUnderTest = new Subject.Builder()
                .principals(new SimplePrincipalCollection(ClothoRealm.ANONYMOUS_USER, "clotho"))
                .authenticated(true)
                .buildSubject();
        setSubject(subjectUnderTest);        
    }
}

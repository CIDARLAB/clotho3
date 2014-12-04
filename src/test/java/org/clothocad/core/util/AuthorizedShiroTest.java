/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import org.apache.shiro.subject.Subject;
import org.clothocad.core.security.AbstractShiroTest;
import org.clothocad.core.security.ServerSubject;
import org.junit.Before;

/**
 *
 * @author spaige
 */
public class AuthorizedShiroTest extends AbstractShiroTest {

    public AuthorizedShiroTest() {
        super();
        Subject subjectUnderTest = new ServerSubject();
        setSubject(subjectUnderTest);     
    }
}

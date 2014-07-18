/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import com.google.inject.Guice;
import com.google.inject.Injector;
import org.apache.shiro.SecurityUtils;
import org.apache.shiro.UnavailableSecurityManagerException;
import org.apache.shiro.mgt.SecurityManager;
import org.apache.shiro.subject.Subject;
import org.apache.shiro.subject.support.SubjectThreadState;
import org.apache.shiro.util.LifecycleUtils;
import org.apache.shiro.util.ThreadState;
import org.clothocad.core.testers.ClothoTestModule;
import org.clothocad.core.util.JongoTestModule;
import org.junit.AfterClass;
 
/**
 * Abstract test case enabling Shiro in test environments.
 */
public abstract class AbstractShiroTest {
 
    private static ThreadState subjectThreadState;
    protected Injector injector;
    
    public AbstractShiroTest() {
        injector = Guice.createInjector(new ClothoTestModule(), new JongoTestModule());
        setSecurityManager(injector.getInstance(SecurityManager.class));
    }
    
   /* @BeforeClass
    public static void beforeClass() {
        //0.  Build and set the SecurityManager used to build Subject instances used in your tests
        //    This typically only needs to be done once per class if your shiro.ini doesn't change,
        //    otherwise, you'll need to do this logic in each test that is different
        //Factory<SecurityManager> factory = new IniSecurityManagerFactory("classpath:test.shiro.ini");

    }*/
    
    //@Before 
    public void beforeTest(){
        //1.  Build the Subject instance for the test to run:
        Subject subjectUnderTest = new Subject.Builder(getSecurityManager()).buildSubject();
 
        //2. Bind the subject to the current thread:
        setSubject(subjectUnderTest);        
    }
    
    //@After
    public void afterTest(){
        clearSubject();
    }
    /**
     * Allows subclasses to set the currently executing {@link Subject} instance.
     *
     * @param subject the Subject instance
     */
    protected void setSubject(Subject subject) {
        clearSubject();
        subjectThreadState = createThreadState(subject);
        subjectThreadState.bind();
    }
 
    protected Subject getSubject() {
        return SecurityUtils.getSubject();
    }
 
    protected ThreadState createThreadState(Subject subject) {
        return new SubjectThreadState(subject);
    }
 
    /**
     * Clears Shiro's thread state, ensuring the thread remains clean for future test execution.
     */
    protected void clearSubject() {
        doClearSubject();
    }
 
    private static void doClearSubject() {
        if (subjectThreadState != null) {
            subjectThreadState.clear();
            subjectThreadState = null;
        }
    }
 
    protected static void setSecurityManager(SecurityManager securityManager) {
        SecurityUtils.setSecurityManager(securityManager);
    }
 
    protected static SecurityManager getSecurityManager() {
        return SecurityUtils.getSecurityManager();
    }
 
    @AfterClass
    public static void tearDownShiro() {
        doClearSubject();
        try {
            SecurityManager securityManager = getSecurityManager();
            LifecycleUtils.destroy(securityManager);
        } catch (UnavailableSecurityManagerException e) {
            //we don't care about this when cleaning up the test environment
            //(for example, maybe the subclass is a unit test and it didn't
            // need a SecurityManager instance because it was using only
            // mock Subject instances)
        }
        setSecurityManager(null);
    }
}

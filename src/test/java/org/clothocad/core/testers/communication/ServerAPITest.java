/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.testers.communication;

import com.google.inject.Guice;
import com.google.inject.Injector;
import org.clothocad.core.layers.communication.ServerSideAPI;
import org.clothocad.core.layers.communication.mind.Mind;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.persistence.mongodb.MongoDBModule;
import org.clothocad.core.testers.ClothoTestModule;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.Ignore;

/**
 *
 * @author spaige
 */
public class ServerAPITest {
    
    private static ServerSideAPI api;
    private static ServerSideAPI unprivilegedUser;
    private static Persistor persistor; 
   
    public ServerAPITest() {}
    
    @BeforeClass
    public static void setUpClass() {
        Injector injector = Guice.createInjector(new ClothoTestModule(), new MongoDBModule());
        persistor = injector.getInstance(Persistor.class);
        Mind mind = new Mind();
        api = new ServerSideAPI(mind,persistor);
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
                persistor.deleteAll();
        //set up test data database
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void get(){
        
    }
    
    @Test
    public void getNonExistent(){
        
    }
    
    @Test
    @Ignore("authz not ready")
    public void getPrivate(){
        
    }
    
    @Test 
    public void set(){
        
    }
    
    @Test
    public void setNonExistent(){
        
    }
    
    @Test
    public void setInvalid(){
        
    }
    
    @Test
    @Ignore("authz not ready")
    public void setPrivate(){
        
    }
    
    @Test 
    public void create(){
        
    }
    
    @Test
    public void createExisting(){
        
    }
    
    @Test 
    @Ignore("authz not ready")
    public  void createWithoutPrivs(){
        
    }
    
    @Test
    public void destroy(){
        
    }
    
    
    @Test
    @Ignore("authz not ready")
    public void destroyWithoutPrivs(){
        
    }
    
    @Test
    public void query(){
        //filter out unseen results
    }
    
    @Test
    public void run(){
        
    }
    
    
    @Test
    public void runNonExistent(){
        
    }
    
    @Test
    public void runWrongArguments(){
        
    }
    
    @Test
    public void runExecutionError(){
        
    }
    
    @Test
    @Ignore("authz not ready")
    public void runWithoutPrivs(){
        
    }
    
}
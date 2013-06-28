/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.testers.communication;

import org.clothocad.core.layers.communication.ServerSideAPI;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author spaige
 */
public class ServerAPITest {
    
    private static ServerSideAPI api;
    private static ServerSideAPI unprivilegedUser;
    
    public ServerAPITest() {
        
    }
    
    @BeforeClass
    public static void setUpClass() {   
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public static void get(){
        
    }
    
    
    @Test
    public static void getNonExistent(){
        
    }
    
    @Test
    public static void getPrivate(){
        
    }
    
    @Test 
    public static void set(){
        
    }
    
    @Test
    public static void setNonExistent(){
        
    }
    
    @Test
    public static void setInvalid(){
        
    }
    
    @Test
    public static void setPrivate(){
        
    }
    
    @Test 
    public static void create(){
        
    }
    
    @Test
    public static void createExisting(){
        
    }
    
    @Test 
    public static void createWithoutPrivs(){
        
    }
    
    @Test
    public static void destroy(){
        
    }
    
    
    @Test
    public static void destroyWithoutPrivs(){
        
    }
    
    @Test
    public static void query(){
        //filter out unseen results
    }
}
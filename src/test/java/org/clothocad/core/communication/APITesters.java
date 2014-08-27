/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.clothocad.core.communication;

import com.google.inject.Guice;
import com.google.inject.Injector;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.apache.shiro.SecurityUtils;
import static org.clothocad.core.ReservedFieldNames.ID;
import org.clothocad.core.datums.Argument;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.execution.ConverterFunction;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.persistence.jongo.JongoModule;
import org.clothocad.core.schema.Schema;
import org.clothocad.core.security.ClothoRealm;
import org.clothocad.core.testers.ClothoTestModule;
import org.clothocad.core.util.JSON;
import org.clothocad.core.util.TestUtils;
import org.clothocad.model.NucSeq;
import org.clothocad.model.Person;
import org.clothocad.model.SimpleSequence;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;

/**
 *
 * @author prashantvaidyanathan
 */
public class APITesters {
    
    private static Router router;
    private static Injector injector;
    private static List<ObjectId> ids;
    private static Persistor persistor;
    
    public APITesters() {
    }

    @BeforeClass
    public static void setUpClass() {
        injector = Guice.createInjector(new ClothoTestModule(), new JongoModule());
        persistor = injector.getInstance(Persistor.class);
        router = injector.getInstance(Router.class);
        org.apache.shiro.mgt.SecurityManager securityManager = injector.getInstance(org.apache.shiro.mgt.SecurityManager.class);
        SecurityUtils.setSecurityManager(securityManager);
        ClothoRealm realm = injector.getInstance(ClothoRealm.class);
        TestUtils.setupTestUsers(realm);
    }

    @AfterClass
    public static void tearDownClass() {
    }

    @Before
    public void setUp() {
        //persistor.deleteAll();
        //ids = TestUtils.setupTestData(persistor);
        injector.getInstance(Persistor.class).deleteAll();
        ids = TestUtils.setupTestData(injector.getInstance(Persistor.class));
        
    }

    @After
    public void tearDown() {
        injector.getInstance(Persistor.class).deleteAll();
        ids = TestUtils.setupTestData(injector.getInstance(Persistor.class));
        
    }
  

    
    private void sendMessage(Message message, ClientConnection connection) throws IOException {
        String stringMessage = JSON.serializeForExternal(message);
        message = JSON.mapper.readValue(stringMessage, Message.class);
        router.receiveMessage(connection, message);
    }
    
    
    public void logintest() throws IOException
    {
        SimpleSequence s = new SimpleSequence("Simple Seq","atgc");
        NucSeq ns = new NucSeq("ATTGGCCTTAAAA");
        //System.out.println("SimpleSeq Class : "+s.getClass().getCanonicalName());
        Argument arguments[];
        Argument arg1 = new Argument("",s.getClass());
        //Argument arg2 = new Argument();
        
        
        ObjectId id1 = new ObjectId();
        String name = "testConv";
        
        TestConnection connection = new TestConnection("LoginTest");
        
        Map<String,Object> simpleseqData = new HashMap<>();
        simpleseqData.put("name", "SimplestSequence");
        simpleseqData.put("schema","org.clothocad.model.SimpleSequence");
        simpleseqData.put("sequence","ATTGGCCTTAAACCC");
        persistor.save(ns);
        
        Map<String,Object> convertparams = new HashMap<>();
       
        //convertparams.put("convertTo", persistor.get(Schema.class, new ObjectId("org.clothocad.model.SimpleSequence")));
         //credentials.put("username", "testuser");
        //credentials.put("password", "password");
         /*final Message message = new Message(
            Channel.createAll,
            new Map[] {simpleseqData},
            "1",
            null
        );*/
       //sendMessage(message, connection);
        Map<String,String> credentials = new HashMap<>();
        credentials.put("username", "maxbates");
        credentials.put("password", "password2");
       
       
        
        sendMessage(new Message(Channel.login, credentials, "1"), connection);
        
        //System.out.println("Created");
        Map<String, Object> operson = persistor.getAsJSON(new ObjectId("clotho.developer.maxbates"));
        System.out.println(operson);
       
    }
    public void createnewusertest() throws IOException
    {
        
        
        TestConnection connection = new TestConnection("LoginTest");
        
        
        Map<String,String> credentials = new HashMap<>();
        credentials.put("username", "vprashant1@live.com");
        credentials.put("password", "mypassword");
        credentials.put("displayname", "PrashantVaidyanathan");
         
        Map<String,String> credentials2 = new HashMap<>();
        credentials2.put("username", "vprashant1@live.com");
        credentials2.put("password", "mypassword");
        credentials2.put("displayname", "PrashantVaidyanathan");
       
        
        Map<String,String> logincred = new HashMap<>();
        logincred.put("username", "vprashant1@gmail.com");
        logincred.put("password", "mypassword");
        
        
        sendMessage(new Message(Channel.createUser, credentials, "1"), connection);
        sendMessage(new Message(Channel.createUser, credentials2, "2"), connection);
        
        sendMessage(new Message(Channel.login, logincred, "3"), connection);
        
        //System.out.println("Created");
        /*Map<String, Object> operson = persistor.getAsJSON(new ObjectId("clotho.developer.maxbates"));
        System.out.println(operson);
        */
    }
    
    
    
    public static void main(String args[]) throws IOException
    {
        APITesters x = new APITesters();
        
        setUpClass();
        x.setUp();
        //x.logintest();
        x.createnewusertest();
    }
    
    
    
}

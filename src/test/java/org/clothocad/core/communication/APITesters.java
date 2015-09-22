/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.clothocad.core.communication;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.clothocad.core.datums.Argument;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.util.JSON;
import org.clothocad.model.Person;
import org.clothocad.model.Sequence;

/**
 *
 * @author prashantvaidyanathan
 */
public class APITesters extends AbstractServerAPITest {

    public APITesters() {
        super();
    }
    
    private void sendMessage(Message message, ClientConnection connection) throws IOException {
        String stringMessage = JSON.serializeForExternal(message);
        message = JSON.mapper.readValue(stringMessage, Message.class);
        router.receiveMessage(connection, message);
    }
    
    
    public void logintest() throws IOException
    {
        Person author = new Person("Temp");
    	Sequence s1 = new Sequence("seq1", "atgc", author);
    	Sequence s2 = new Sequence("seq2","ATTGGCCTTAAAA", author);
        //System.out.println("SimpleSeq Class : "+s.getClass().getCanonicalName());
        Argument arguments[];
        Argument arg1 = new Argument("",s1.getClass());
        //Argument arg2 = new Argument();
        
        
        ObjectId id1 = new ObjectId();
        String name = "testConv";
        
        TestConnection connection = new TestConnection("LoginTest");
        
        Map<String,Object> simpleseqData = new HashMap<>();
        simpleseqData.put("name", "SimplestSequence");
        simpleseqData.put("schema","org.clothocad.model.Sequence");
        simpleseqData.put("sequence","ATTGGCCTTAAACCC");
        persistor.save(s2);
        
        Map<String,Object> convertparams = new HashMap<>();
       
        //convertparams.put("convertTo", persistor.get(Schema.class, new ObjectId("org.clothocad.model.Sequence")));
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
        
        Map<String,String> logincred2 = new HashMap<>();
        logincred2.put("username", "vprashant1@live.com");
        logincred2.put("password", "anotherpass");
        
         Map<String,String> logincred3 = new HashMap<>();
        logincred3.put("username", "vprashant1@live.com");
        logincred3.put("password", "mypassword");
        
        
        sendMessage(new Message(Channel.createUser, credentials, "1"), connection);
        sendMessage(new Message(Channel.createUser, credentials2, "2"), connection);
        
        sendMessage(new Message(Channel.login, logincred, "3"), connection);
        sendMessage(new Message(Channel.updatePassword, logincred2, "4"), connection);
        sendMessage(new Message(Channel.login, logincred3, "5"), connection);
        
        //System.out.println("Created");
        /*Map<String, Object> operson = persistor.getAsJSON(new ObjectId("clotho.developer.maxbates"));
        System.out.println(operson);
        */
    }
}

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.clothocad.core.execution;

import com.fasterxml.jackson.core.JsonParseException;
import com.google.inject.Guice;
import com.google.inject.Injector;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import static org.clothocad.core.ReservedFieldNames.ID;
import org.clothocad.core.communication.Channel;
import org.clothocad.core.communication.ClientConnection;
import org.clothocad.core.communication.Message;
import org.clothocad.core.communication.Router;
import org.clothocad.core.communication.ServerSideAPI;
import org.clothocad.core.communication.TestConnection;
import org.clothocad.core.datums.Argument;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.persistence.jongo.JongoModule;
import org.clothocad.core.schema.Schema;
import org.clothocad.core.testers.ClothoTestModule;
import org.clothocad.core.util.JSON;
import org.clothocad.core.util.TestUtils;
import org.clothocad.model.NucSeq;
import org.clothocad.model.SimpleSequence;
import org.junit.After;
import org.junit.AfterClass;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

/**
 *
 * @author prashantvaidyanathan
 */
public class ConverterTest {
     private static ServerSideAPI api;
    private static ServerSideAPI unprivilegedUser;
    private static Persistor persistor;
    private static List<ObjectId> ids;
    private static Mind mind;
    private static Router router;

    public ConverterTest() {
    }

    @BeforeClass
    public static void setUpClass() {
       Injector injector = Guice.createInjector(new ClothoTestModule(), new JongoModule());
        persistor = injector.getInstance(Persistor.class);
        router = injector.getInstance(Router.class);
        mind = new Mind();
        api = new ServerSideAPI(mind, persistor, router, null);
        persistor.connect();
        mind.setConnection(new TestConnection("test"));
    }

    @AfterClass
    public static void tearDownClass() {
    }

    @Before
    public void setUp() {
        persistor.deleteAll();
        ids = TestUtils.setupTestData(persistor);
    }

    @After
    public void tearDown() {
    }
    // TODO add test methods here.
    // The methods must be annotated with annotation @Test. For example:
    //
    // @Test
    // public void hello() {}

    
    
    

    private void sendMessage(Message message, ClientConnection connection) throws IOException {
        String stringMessage = JSON.serializeForExternal(message);
        message = JSON.mapper.readValue(stringMessage, Message.class);
        router.receiveMessage(connection, message);
    }
    
    public void converttest() throws IOException
    {
        SimpleSequence s = new SimpleSequence("Simple Seq","atgc");
        NucSeq ns = new NucSeq("ATTGGCCTTAAAA");
        System.out.println("SimpleSeq Class : "+s.getClass().getCanonicalName());
        Argument arguments[];
        Argument arg1 = new Argument("",s.getClass());
        //Argument arg2 = new Argument();
        
        ConverterFunction convfunc1 = new ConverterFunction();
        ConverterFunction convfunc2 = new ConverterFunction();
        
        convfunc1.convertTo = persistor.get(Schema.class, new ObjectId("org.clothocad.model.NucSeq"));
        convfunc1.convertFrom = persistor.get(Schema.class, new ObjectId("org.clothocad.model.SimpleSequence"));
        convfunc1.setName("Func1");
        convfunc2.convertTo = persistor.get(Schema.class, new ObjectId("org.clothocad.model.SimpleSequence"));
        convfunc2.convertFrom = persistor.get(Schema.class, new ObjectId("org.clothocad.model.NucSeq"));
        convfunc2.setName("Func2");
        Function func1 = new Function();
        Function func2 = new Function();
        
        convfunc1.convFunction = func1;
        convfunc2.convFunction = func2;
        
        persistor.save(convfunc1);
        persistor.save(convfunc2);
        
        
       
        ObjectId id1 = new ObjectId();
        ObjectId id2 = new ObjectId();
        String name = "testConv";
        Function convFunc = new Function();
        convFunc.setCode("function(sequence) { return sequence.split('').reverse().join('');};");
        convFunc.setLanguage(Language.JAVASCRIPT);
        convFunc.setName("testConv");
        
        System.out.println("Hello!");
        TestConnection connection = new TestConnection("ConverterTest");
        
        /*
        Map<String,Object> simpleseqData = new HashMap<>();
        simpleseqData.put("name", "SimplestSequence");
        simpleseqData.put("schema","org.clothocad.model.SimpleSequence");
        simpleseqData.put("sequence","ATTGGCCTTAAACCC");
        
        Map<String,Object> convertparams = new HashMap<>();
        convertparams.put("name","convertAct");
        convertparams.put("convert", simpleseqData);
        convertparams.put("convertTo", "org.clothocad.model.SimpleSequence");
        
        */
        
        /*final Message message = new Message(
            Channel.createAll,
            new Map[] {simpleseqData},
            "1",
            null
        );*/
        //sendMessage(message, connection);
        //sendMessage(new Message(Channel.convert, convertparams, "2"), connection);
        
       
        
        Map<String,Object> convfunction1 = new HashMap<>();
        convfunction1.put("schema", "org.clothocad.core.execution.ConverterFunction");
        convfunction1.put("name","NucSeqtoSimpleSeq");
        //convfunction1.put("convFunction", convFunc);
        convfunction1.put("convertFrom", "org.clothocad.model.NucSeq");
        convfunction1.put("convertTo", "org.clothocad.mode.SimpleSequence");
        convfunction1.put(ID,id1.toString());
        
        
        Map<String,Object> convfunction2 = new HashMap<>();
        convfunction2.put("schema", "org.clothocad.core.execution.ConverterFunction");
        convfunction2.put("name","SimpleSeqtoNucSeq");
        convfunction2.put("convertFrom", "org.clothocad.mode.SimpleSequence");
        convfunction2.put("convertTo", "org.clothocad.model.NucSeq");
        convfunction2.put(ID,id2.toString());
        
        //ObjectId createdId1 = api.create(convfunction1);
        //ObjectId createdId2 = api.create(convfunction2);
        System.out.println("Created");
        
        
        
        Map<String, Object> query = new HashMap<>();
        query.put("schema", "org.clothocad.core.execution.ConverterFunction");
        
        List<Map<String, Object>> results = api.query(query);
        Set<String> names = new HashSet();
        
        Collection<ConverterFunction> convlist = persistor.getAll(ConverterFunction.class);
        for(ConverterFunction xconvfunc : convlist)
        {
            if(xconvfunc.convertTo.equals(persistor.get(Schema.class, new ObjectId("org.clothocad.model.NucSeq"))))
            {
                if(xconvfunc.convertFrom.equals(persistor.get(Schema.class, new ObjectId("org.clothocad.model.SimpleSequence"))))
                {
                    System.out.println("We found the converter function you were looking for");
                    System.out.println("Name of this function" + xconvfunc.getName());
                }
            }
        }
        /*
        for (Map<String, Object> result : results) {
            Object convF = result.get("convertFrom");
            Object convT = result.get("convertTo");
            if(convF != null && convT != null)
            {
                System.out.println("Not null");
                
                if(convF.equals(persistor.get(Schema.class, new ObjectId("org.clothocad.model.NucSeq"))))
                {
                    if(convF.equals(persistor.get(Schema.class, new ObjectId("org.clothocad.model.SimpleSequence"))))
                    {
                        System.out.println("Found the converter we were looking for!!");
                        System.out.println("Name of this converter is : "+ result.get("name").toString());
                    }
                }
            }
            
            //names.add(result.get("name").toString());
            System.out.println(result.get("name").toString());
        }
        */
        
        /*final Message message = new Message(
            Channel.createAll,
            new Map[] {convfunction1,convfunction2},
            "2",
            null
        );
        sendMessage(message, connection);
        System.out.println("Converter Created!");
        */
       
        /*
        Map querystuff = new HashMap<>();
        //querystuff.put("tokens", new ArrayList());
        querystuff.put("schema","org.clothocad.core.execution.ConverterFunction");
        //querystuff.put("query", "clotho.get(\"org.clothocad.core.execution.ConverterFunction\")");
        sendMessage(new Message(Channel.query, querystuff, "3"), connection);
        Map<String,Object> data = (Map) connection.messageDataByChannelAndId.get(Channel.submit.name()+"3");
        for(Map.Entry<String,Object> item : data.entrySet())
        {
            System.out.println("Key : " + item.getKey() + ":"+item.getValue());
        }
        */
        //System.out.println(x.convertTo.toString());
    }
    
    public static void main(String args[]) throws IOException
    {
        ConverterTest x = new ConverterTest();
        
        setUpClass();
        x.setUp();
        x.converttest();
    }
}

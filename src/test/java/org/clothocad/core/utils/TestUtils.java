/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.utils;

import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.MappingIterator;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.ObjectReader;
import com.google.inject.Guice;
import com.google.inject.Injector;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import lombok.extern.slf4j.Slf4j;
import org.bson.types.ObjectId;
import org.clothocad.core.layers.communication.ServerSideAPI;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.persistence.mongodb.MongoDBModule;
import org.clothocad.core.security.ClothoRealm;
import org.clothocad.core.security.CredentialStore;
import org.clothocad.core.testers.ClothoTestModule;
import org.clothocad.model.FreeForm;
import org.clothocad.model.Institution;
import org.clothocad.model.Lab;
import org.clothocad.model.Part;
import org.clothocad.model.Part.PartFunction;
import org.clothocad.model.Person;

/**
 *
 * @author spaige
 */
@Slf4j
public class TestUtils {

    private static void importTestJSON(Persistor persistor) {
        ServerSideAPI api = new DummyAPI(persistor);
        ObjectReader reader = new ObjectMapper().reader(Map.class);
        for (File child : new File("src/test/resources/testData").listFiles()){
            try {
                MappingIterator<Map> it = reader.readValues(child);
                while (it.hasNext()){
                    api.create(it.next());
                }
            } catch (JsonProcessingException ex) {
                log.warn("Could not process {} as JSON", child.getAbsolutePath());
            } catch (IOException ex) {
                log.warn("Could not open {}", child.getAbsolutePath());
            }
        }
    }
    //not sure if this should be static 
    
    private Injector injector;

    public <T> T getA(Class<T> type) {
        if (injector == null){
            injector = getDefaultTestInjector();
        }
        return injector.getInstance(type);
    }
    
    public static Injector getDefaultTestInjector(){
        return Guice.createInjector(new ClothoTestModule(), new MongoDBModule());
    }
    
    public static List<ObjectId> setupTestData(Persistor persistor){
        importTestJSON(persistor);
        
        Institution i = new Institution("Test institution", "Townsville", "Massachusetts", "United States of America");
        Lab lab = new Lab(i, null, "Test Lab", "College of Testing", "8 West Testerfield");
        Person person = new Person("Test Person", lab, null);
        lab.setPI(person);
        
        Part part1 = Part.generateBasic("Test Part 1", "the first test part", "AAAAAAAAAAAAAAAAAAA", new FreeForm(), person);
        part1.setType(PartFunction.CDS);
        
        Part part2 = Part.generateBasic("Test Part 2", "the second test part", "TTTTTTTTTTTTTTTTTT", new FreeForm(), person);
        Part part3 = Part.generateComposite(Arrays.asList(part1, part2), new Object(), new FreeForm(), person, "Test Part 3", "parts 1 and 2 jammed together");
                
        Map<String,Object> eugenePart = new HashMap();
        eugenePart.put("Name", "B0015");
        eugenePart.put("schema", "eugene.dom.components.Part");
        eugenePart.put("PartType", "Terminator");
        eugenePart.put("Sequence", "CCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGGTGAACGCTCTCTACTAGAGTCACACTGGCTCACCTTCGGGTGGGCCTTTCTGCGTTTATA");
        eugenePart.put("Pigeon", "t B0015");
        
        persistor.save(part1);
        persistor.save(part2);
        persistor.save(part3);
        ObjectId eugeneID = persistor.save(eugenePart);
        
        return Arrays.asList(part1.getUUID(), part2.getUUID(), part3.getUUID(), eugeneID);
    }
    
    public static void setupTestUsers(ClothoRealm realm){
        realm.addAccount("testuser", "password");
    }
}

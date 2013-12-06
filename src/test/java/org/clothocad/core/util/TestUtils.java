/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.MappingIterator;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.ObjectReader;
import com.google.inject.Guice;
import com.google.inject.Injector;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.persistence.EntityNotFoundException;
import lombok.extern.slf4j.Slf4j;
import org.bson.types.ObjectId;
import org.clothocad.core.communication.ServerSideAPI;
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

    public static void importTestJSON(Persistor persistor) {
        importTestJSON("src/test/resources/testData", persistor, false);
    }

    public static void importTestJSON(String path, Persistor persistor, boolean overwrite) {
        ServerSideAPI api = new DummyAPI(persistor);
        ObjectReader reader = new ObjectMapper().reader(Map.class);
        
        List<Map> objects = new ArrayList<>();
        for (File child : new File(path).listFiles()) {
            if (!child.getName().endsWith(".json")) {
                continue;
            }
            try {
                MappingIterator<Map> it = reader.readValues(child);
                while (it.hasNext()) {
                    objects.add(it.next());
                }
                
            } catch (JsonProcessingException ex) {
                log.warn("Could not process {} as JSON", child.getAbsolutePath());
            } catch (IOException ex) {
                log.warn("Could not open {}", child.getAbsolutePath());
            }
        }
        
        while (objects.size() > 0) {
            int prevSize = objects.size();
            List<Map> newObjects = new ArrayList();

            for (Map obj : objects) {
                if (!overwrite) {
                    try {
                        persistor.resolveSelector(obj.get("name").toString(), false);
                        continue;

                    } catch (EntityNotFoundException e) {
                    }
                }
                try {
                    ObjectId result = api.create(obj);

                    if (overwrite && result == null) {
                        api.set(obj);
                    }
                } catch (RuntimeException e) {
                    newObjects.add(obj);
                } catch (ClassCircularityError e){
                    e.getMessage();
                }
            }

            objects = newObjects;
            if (objects.size() >= prevSize) {
                log.error("Could not load some files: {}", objects.toString());
                return;
            }

        }
    }
    //not sure if this should be static 
    private Injector injector;

    public <T> T getA(Class<T> type) {
        if (injector == null) {
            injector = getDefaultTestInjector();
        }
        return injector.getInstance(type);
    }

    public static Injector getDefaultTestInjector() {
        return Guice.createInjector(new ClothoTestModule(), new MongoDBModule());
    }

    public static List<ObjectId> setupTestData(Persistor persistor) {
        importTestJSON(persistor);

        Institution i = new Institution("Test institution", "Townsville", "Massachusetts", "United States of America");
        Lab lab = new Lab(i, null, "Test Lab", "College of Testing", "8 West Testerfield");
        Person person = new Person("Test Person",null);
        lab.setPI(person);

        Part part1 = Part.generateBasic("Test Part 1", "the first test part", "AAAAAAAAAAAAAAAAAAA", new FreeForm(), person);
        part1.setType(PartFunction.CDS);

        Part part2 = Part.generateBasic("Test Part 2", "the second test part", "TTTTTTTTTTTTTTTTTT", new FreeForm(), person);
        Part part3 = Part.generateComposite(Arrays.asList(part1, part2), new Object(), new FreeForm(), person, "Test Part 3", "parts 1 and 2 jammed together");

        Map<String, Object> eugenePart = new HashMap();
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

    public static void setupTestUsers(ClothoRealm realm) {
        realm.addAccount("testuser", "password");
    }
}

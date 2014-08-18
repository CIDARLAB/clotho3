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
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.persistence.EntityNotFoundException;
import lombok.extern.slf4j.Slf4j;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.persistence.jongo.JongoModule;
import org.clothocad.core.security.ClothoRealm;
import org.clothocad.core.testers.ClothoTestModule;
import org.clothocad.model.FreeForm;
import org.clothocad.model.Institution;
import org.clothocad.model.Lab;
import org.clothocad.model.Part;
import org.clothocad.model.Part.PartFunction;
import org.clothocad.model.Person;
import static org.clothocad.core.ReservedFieldNames.*;
import org.clothocad.core.persistence.ClothoConnection;

/**
 *
 * @author spaige
 */
@Slf4j
public class TestUtils {

    public static void importTestJSON(Persistor persistor) {
        importTestJSON(Paths.get("src","test","resources","testData"), persistor.getConnection(), true);
    }

    public static void importTestJSON(Path path, ClothoConnection connection, boolean overwrite) {
        ObjectReader reader = new ObjectMapper().reader(Map.class);
        
        List<Map> objects = new ArrayList<>();
        try {
            for (Path child : Files.newDirectoryStream(path)) {
                if (!child.toString().endsWith(".json")) {
                    continue;
                }
                try {
                    MappingIterator<Map> it = reader.readValues(child.toFile());
                    while (it.hasNext()) {
                        objects.add(it.next());
                    }
                    
                } catch (JsonProcessingException ex) {
                    log.warn("Could not process {} as JSON", child.toAbsolutePath());
                } catch (IOException ex) {
                    log.warn("Could not open {}", child.toAbsolutePath());
                }
            }
        } catch (IOException ex) {
            log.warn("Could not open {}", path.toAbsolutePath());
            throw new RuntimeException(ex);
        }
        for (Map obj : objects) {
            if (!overwrite) {
                try {
                    if (obj.containsKey(ID) && connection.exists(new ObjectId(obj.get(ID)))) {
                        //don't overwrite things
                        continue;
                    }

                } catch (EntityNotFoundException e) {
                }
            }
            connection.save(obj);
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
        return Guice.createInjector(new ClothoTestModule(), new JongoModule());
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
        
        Function dummyPackager = DummyPackager.createDummyPackager();
        persistor.save(dummyPackager);
        
        ObjectId eugeneID = persistor.save(eugenePart);

        return Arrays.asList(part1.getId(), part2.getId(), part3.getId(), eugeneID);
    }

    public static void setupTestUsers(ClothoRealm realm) {
        realm.addAccount("testuser", "password");
    }
    
    public static Map<String, Object> serializeForExternalAsMap(Object o){
        String serialized = JSON.serializeForExternal(o);
        try {
            return JSON.deserializeObjectToMap(serialized);
        } catch (IOException ex) {
            ex.printStackTrace();
            return null;
        }
    }
}

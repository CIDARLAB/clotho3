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
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.concurrent.Callable;
import javax.persistence.EntityNotFoundException;
import lombok.extern.slf4j.Slf4j;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.Persistor;
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
import org.clothocad.core.security.SecurityModule;

/**
 *
 * @author spaige
 */
@Slf4j
public class TestUtils {

    public static void importTestJSON(Persistor persistor) {
        importJSONFromDirectory(Paths.get("src","test","resources","testData"), persistor.getConnection(), null, true, false);
    }

    public static void importJSONFromDirectory(Path path, ClothoConnection connection, ClothoRealm realm, boolean overwrite, boolean allPublic) {        
        List<Map<String,Object>> objects;
        try {
            objects = readJSON(Files.newDirectoryStream(path));
        } catch (IOException ex) {
            log.error("Could not open {}", path.toAbsolutePath());
            throw new RuntimeException(ex);
        }
        saveJSONObjects(objects, connection, realm, overwrite, allPublic);
    }
    
    public static List<Map<String,Object>> readJSON(Iterable<Path> paths){
        ObjectReader reader = new ObjectMapper().reader(Map.class);
        List<Map<String,Object>> output = new ArrayList<>();
        for (Path child : paths) {
            if (!child.toString().endsWith(".json")) {
                continue;
            }
            try {
                MappingIterator<Map> it = reader.readValues(child.toFile());
                while (it.hasNext()) {
                    output.add(it.next());
                }

            } catch (JsonProcessingException ex) {
                log.warn("Could not process {} as JSON", child.toAbsolutePath());
            } catch (IOException ex) {
                log.warn("Could not open {}", child.toAbsolutePath());
            }
        }
        return output;
    }
    
    public static void saveJSONObjects(Iterable<Map<String, Object>> objects,
            ClothoConnection connection, ClothoRealm realm, boolean overwrite, boolean allPublic) {
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
            if (allPublic) {
                realm.setPublic(new ObjectId(obj.get(ID)));
            }
        }
    }
    
    public static void importListedJSON(Path targetDirectory, Path idsListFile,
            ClothoConnection connection, ClothoRealm realm){
        
        List<String> ids;
        
        List<Path> paths = new ArrayList<>();
        try {
            Scanner scanIds = new Scanner(idsListFile);
            while (scanIds.hasNextLine()){
                paths.add(targetDirectory.resolve(scanIds.nextLine()+".json"));
            }
        } catch (IOException ex) {
            log.error("Could not open {}", idsListFile.toAbsolutePath());
        }
        
        List<Map<String,Object>> objects = readJSON(paths);
        saveJSONObjects(objects, connection, realm, true, true);
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
        return Guice.createInjector(new ClothoTestModule(), new SecurityModule(), new JongoTestModule());
    }
    
    public static void setupAuthoringTestData(Persistor persistor)
    {
        importTestJSON(persistor);
        Person newperson = new Person("testuser");
        Map<String, Object> newperson2 = new HashMap();
        newperson2.put("displayname", "maxbates");
        //newperson2.put("rawPassword", "password2");
        newperson2.put("description", "Maxwell Bates graduated from UC Berkeley in 2012 with a B.S. in Bioengineering and focus on Genetic Engineering. He began working on the Clotho project in the summer prior to his graduation, and has been on board as the primary front end developer since. Prior to Clotho, he worked with Dr. Eric Topol at the Scripps Translational Science Institute. His interests revolve around the empowerment of individuals through online education, and the empowerment of patients toward personalizing medicine.");
        newperson2.put("current", true);
        newperson2.put("currentLocation", "Oakland, CA");
        newperson2.put("dateCreated", "2012-06-01");
        newperson2.put("emailAddress", "maxbates@gmail.com");
        newperson2.put("endDate", "");
        newperson2.put("givenName", "Maxwell");
        newperson2.put("icon", "images/people/MaxwellBates.jpg");
        newperson2.put("imgUrl", "");
        newperson2.put("isManagement", false);
        newperson2.put("lab", "Anderson");
        newperson2.put("name", "maxbates");
        newperson2.put("nickName", "Max");
        newperson2.put("role", "Programmer");
        newperson2.put("schema", "org.clothocad.model.Person");
        newperson2.put("social", "");
        newperson2.put("startDate", "2012-06-01");
        newperson2.put("surName", "Bates");
        newperson2.put("title", "Front End Developer");
        newperson2.put("id", "clotho.testuser.maxbates");
        
        //persistor.save(newperson);
        //persistor.save(newperson2);
        
    }

    public static List<ObjectId> setupTestData(Persistor persistor, ClothoRealm realm) {
        importTestJSON(persistor);
        importListedJSON(Paths.get("src", "test", "resources", "authoredJSON"), Paths.get("src", "test", "resources", "requiredAuthoredResources.txt"), persistor.getConnection(), realm);

        Institution i = new Institution("Test institution", "Townsville", "Massachusetts", "United States of America");
        Lab lab = new Lab(i, null, "Test Lab", "College of Testing", "8 West Testerfield");
        Person person = new Person("Test Person");
        lab.setPI(person);
        
        
        Person newperson = new Person("testuser");
        Map<String, Object> newperson2 = new HashMap();
        newperson2.put("displayname", "maxbates");
        //newperson2.put("rawPassword", "password2");
        newperson2.put("description", "Maxwell Bates graduated from UC Berkeley in 2012 with a B.S. in Bioengineering and focus on Genetic Engineering. He began working on the Clotho project in the summer prior to his graduation, and has been on board as the primary front end developer since. Prior to Clotho, he worked with Dr. Eric Topol at the Scripps Translational Science Institute. His interests revolve around the empowerment of individuals through online education, and the empowerment of patients toward personalizing medicine.");
        newperson2.put("current", true);
        newperson2.put("currentLocation", "Oakland, CA");
        newperson2.put("dateCreated", "2012-06-01");
        newperson2.put("emailAddress", "maxbates@gmail.com");
        newperson2.put("endDate", "");
        newperson2.put("givenName", "Maxwell");
        newperson2.put("icon", "images/people/MaxwellBates.jpg");
        newperson2.put("imgUrl", "");
        newperson2.put("isManagement", false);
        newperson2.put("lab", "Anderson");
        newperson2.put("name", "maxbates");
        newperson2.put("nickName", "Max");
        newperson2.put("role", "Programmer");
        newperson2.put("schema", "org.clothocad.model.Person");
        newperson2.put("social", "");
        newperson2.put("startDate", "2012-06-01");
        newperson2.put("surName", "Bates");
        newperson2.put("title", "Front End Developer");
        newperson2.put("id", "clotho.testuser.maxbates");

        persistor.save(newperson);
        //persistor.save(newperson2);

        
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
        realm.addAccount("maxbates","password2");
        
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
    
    public static class SetupTestData implements Callable {
        public SetupTestData(Persistor persistor, ClothoRealm realm){
            this.persistor = persistor;
            this.realm = realm;
        }
        
        private Persistor persistor;
        private ClothoRealm realm;
        
        @Override
        public Object call() throws Exception {
            setupTestData(persistor, realm);
            return persistor;
        }
        
    }
    
    public static class SetupTestRealm implements Callable {
        public SetupTestRealm(ClothoRealm realm){
            this.realm = realm;
        }
        private ClothoRealm realm;

        @Override
        public Object call() throws Exception {
            setupTestUsers(realm);
            return realm;
        }
    }
}

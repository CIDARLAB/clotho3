/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import com.google.inject.Injector;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.schema.BuiltInSchema;
import org.reflections.Reflections;

/**
 *
 * @author spaige
 */
public class BuiltInWriter {
    private static Persistor persistor;
    private static List<BuiltInSchema> schemas;
    
    public static void main(String[] args) throws IOException {
        
        Injector injector = TestUtils.getDefaultTestInjector();
        persistor = injector.getInstance(Persistor.class);
        
        schemas = new ArrayList<>();
        
        Reflections models = new Reflections("org.clothocad");

        for (Class<? extends ObjBase> c : models.getSubTypesOf(ObjBase.class)){
            if (c.getSuperclass() == ObjBase.class){
                makeBuiltIns(c, null, models);
            }
        }
        
        Path file = Paths.get("builtins.json");
        try (BufferedWriter writer = Files.newBufferedWriter(file, StandardCharsets.UTF_8)){
            writer.write(JSON.serializeForExternal(schemas));
        } catch (IOException ex) {
            ex.printStackTrace();
        }        
    }
    
    private static void makeBuiltIns(Class<? extends ObjBase> c, BuiltInSchema superSchema, Reflections ref){
        BuiltInSchema builtIn = new BuiltInSchema(c, superSchema);
        builtIn.setId(new ObjectId());
        schemas.add(builtIn);
        
        for (Class<? extends ObjBase> subClass : ref.getSubTypesOf(c)){
            if (subClass.getSuperclass() == c){
                makeBuiltIns(subClass, builtIn, ref);                
            }
        }
    }
}

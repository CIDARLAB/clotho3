/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import com.google.inject.Guice;
import com.google.inject.Injector;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import org.bson.types.ObjectId;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.persistence.jongo.JongoModule;
import org.clothocad.core.schema.BuiltInSchema;
import org.clothocad.core.schema.Schema;
import org.clothocad.core.testers.ClothoTestModule;
import org.reflections.Reflections;

/**
 *
 * @author spaige
 */
public class BuiltInWriter {
    private static Persistor persistor;
    private static List<BuiltInSchema> schemas;
    
    public static void main(String[] args){
        
        Injector injector = Guice.createInjector(
                new ClothoTestModule(),
                new JongoModule());
        persistor = injector.getInstance(Persistor.class);
        
        schemas = new ArrayList<>();
        
        List<Map<String,Object>> data = new ArrayList<>();
        Reflections models = new Reflections("org.clothocad");

        for (Class<? extends ObjBase> c : models.getSubTypesOf(ObjBase.class)){
            if (c.getSuperclass() == ObjBase.class){
                makeBuiltIns(c, null, models);
            }
        }
        
        for (BuiltInSchema schema : schemas){
            data.add(persistor.toJSON(schema));
        }
        
        String content = JSON.serialize(data);
        Path file = Paths.get("builtins.json");
        try (BufferedWriter writer = Files.newBufferedWriter(file, StandardCharsets.UTF_8)){
            writer.write(content);
        } catch (IOException ex) {
            ex.printStackTrace();
        }        
    }
    
    private static void makeBuiltIns(Class<? extends ObjBase> c, Schema superSchema, Reflections ref){
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

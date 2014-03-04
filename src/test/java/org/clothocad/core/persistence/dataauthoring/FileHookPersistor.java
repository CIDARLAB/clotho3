/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence.dataauthoring;

import com.google.inject.Singleton;
import com.google.inject.name.Named;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.nio.file.Files;
import java.util.Collection;
import java.util.Map;
import java.util.Set;
import javax.inject.Inject;
import javax.validation.ConstraintViolationException;
import lombok.extern.slf4j.Slf4j;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.ClothoConnection;
import org.clothocad.core.persistence.OverwriteConfirmationException;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.util.JSON;

/**
 *
 * @author spaige
 */
@Singleton
@Slf4j
public class FileHookPersistor extends Persistor{
    
    protected Path storageFolder;
    
    @Inject
    public FileHookPersistor(final ClothoConnection connection, @Named("storagefolder") Path storageFolder) throws IOException{
        super(connection, false);
        this.storageFolder = storageFolder;
        if (!Files.exists(storageFolder)) Files.createDirectories(storageFolder);
    }
    
    @Override
    public void initializeBuiltInSchemas(){
        super.initializeBuiltInSchemas();
    }
    
    public void writeToFile(Collection<ObjBase> objs){
        for (ObjBase obj: objs){
            writeToFile(obj);
        }
    }
    
    public void writeToFile(ObjBase o){
        Path file = storageFolder.resolve(o.getId().toString()+".json");
        if (file == null){
            log.warn("No id provided for json object: \n{}", o);
        }
        try (BufferedWriter writer = Files.newBufferedWriter(file, StandardCharsets.UTF_8)){
            writer.write(JSON.serializeForExternal(o, true));
        } catch (IOException ex) {
            log.error("Error serializing object {}: {}:", o, ex);
        }
    }
    
    public void writeToFile(Map<String,Object> json){
        Path file = storageFolder.resolve(json.get("id").toString()+".json");
        String content = JSON.serializeJSONMapForExternal(json, true);
        if (file == null){
            log.warn("No id provided for json object: \n{}", content);
        }
        try (BufferedWriter writer = Files.newBufferedWriter(file, StandardCharsets.UTF_8)){
            writer.write(content);
        } catch (IOException ex) {
            log.error("Error serializing object {}: {}:", json, ex);
        }
    }
    
    @Override
    public ObjectId save(ObjBase obj, boolean overwrite) {
        ObjectId returnvalue = super.save(obj, overwrite);
        Set<ObjBase> objs = getObjBaseSet(obj);
        writeToFile(objs);
        return returnvalue;
    }

    @Override
    public ObjectId save(Map<String, Object> data) throws ConstraintViolationException, OverwriteConfirmationException {
        String id; 
        
        if (!data.containsKey("id")){
            id = new ObjectId().toString();
            data.put("id", id);
        } 
        if (data.containsKey("id") && data.get("id") instanceof ObjectId){
            id = ((ObjectId) data.get("id")).toString();
            data.put("id", id);
        }
        writeToFile(data);
        return super.save(data); //To change body of generated methods, choose Tools | Templates.
    }
    
    
}

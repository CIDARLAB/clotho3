/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.layers.persistence;

import javax.inject.Inject;
import org.bson.types.ObjectId;
import org.clothocad.core.aspects.Persistor;
import org.clothocad.core.schema.Schema;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 *
 * @author spaige
 * 
 * For now, is conventional classloader subclass
 * 
 * 
 * later, if permgen memory becomes a big issue, try http://stackoverflow.com/questions/88235/dealing-with-java-lang-outofmemoryerror-permgen-space-error?rq=1
 * or possibly anonymousclassloader hackery
 * 
 * ...can I modify addClass to make it a no-op through asm?
 * 
 */
public class DBClassLoader extends ClassLoader {
    final static Logger logger = LoggerFactory.getLogger(DBClassLoader.class);
    
    @Inject
    public DBClassLoader(Persistor p){
        super();
        persistor = p; 
    }
    
    Persistor persistor;
    
    /*@Override 
    public Class<?> loadClass(String uuid) throws ClassNotFoundException{
        //check cache
            //fake override of findLoadedClass
        //check parent
            //copy from classloader
        
        //get bytes from db
        Schema schema = persistor.get(new ObjectId(uuid));
        
        
        //get classdef from anonymousclassloader
        
        //put in cache
        
        //return class
        
    }*/
    
    @Override
    public Class<?> findClass(String id) throws ClassNotFoundException{
        Schema s = persistor.get(Schema.class, new ObjectId(id));
        if (s == null){
            logger.error("Could not find schema with id {} in database.", id);
            throw new ClassNotFoundException();
        }
        byte[] classData = s.getClassData();
        String className = s.getBinaryName();
        logger.debug("Loading schema with id {} named {} as {} ...", id, s.getName(), className);
        return this.defineClass(className, classData, 0, classData.length);
    }
}

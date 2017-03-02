/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.webserver.jetty;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.model.Annotation;
import org.clothocad.model.BasicModule;
import org.clothocad.model.BioDesign;
import org.clothocad.model.Feature;
import org.clothocad.model.Module;
import org.clothocad.model.Part;
import org.clothocad.model.Person;
import org.clothocad.model.Sequence;

/**
 *
 * @author David
 */
public class ConvenienceMethods {
    
    //Really only used in the REST API, so have REST API pass the persistor to the functions.
    
    public static ObjectId createPart(Persistor persistor, String name, String author) {
        Person auth = new Person(author);
        Part part = new Part(name, auth);
        
        BioDesign design = new BioDesign(name, auth);
        design.addPart(part);
        
        ObjectId id = persistor.save(design);
        
        return id;
        
    }


    public static ObjectId createPart(Persistor persistor, String name, Map<String, String> parameters, String author) {
        if (parameters.isEmpty()) {
            return createPart(persistor, name, author);
        }
        boolean bRole = false, bSeq = false;
        String role = "", sequence = "";

        for (String field : parameters.keySet()) {
            if (field.equalsIgnoreCase("role")) {
                bRole = true;
                role = parameters.get(field);
            }

            if (field.equalsIgnoreCase("sequence")) {
                bSeq = true;
                sequence = parameters.get(field);
            }
        }
        
        Person auth = new Person(author);

        if (bSeq && bRole) {

            Sequence seq = new Sequence(name, sequence, auth);

            Part part = new Part(name, seq, auth);

            Feature feat = new Feature(name, role, auth);
            feat.setSequence(seq);

            Annotation annotation = seq.createAnnotation(name, 1, sequence.length(), true, auth);
            annotation.setFeature(feat);

            Set<Feature> feats = new HashSet<>();
            feats.add(feat);

            Module bMod = new BasicModule(name, role, feats, auth);

            BioDesign design = new BioDesign(name, auth);
            design.addPart(part);
            design.setModule(bMod);
            ObjectId id = persistor.save(design);

            return id;
        }
        else if(bSeq)
        {
            Sequence seq = new Sequence(name, sequence, auth);
            seq.createAnnotation(name, 1, sequence.length(), true, auth);

            Part part = new Part(name, seq, auth);

            BioDesign design = new BioDesign(name, auth);
            design.addPart(part);

            ObjectId id = persistor.save(design);
            
            return id;
            
        }
        else if(bRole)
        {
            Part part = new Part(name, auth);
            Feature feat = new Feature(name, role, auth);
            Annotation annotation = new Annotation();
            annotation.setFeature(feat);
            annotation.setAuthor(auth);
            annotation.setName(name);
            annotation.setForwardStrand(true);
            annotation.setStart(0);
            annotation.setEnd(0);

            Set<Feature> feats = new HashSet<>();
            feats.add(feat);

            Module bMod = new BasicModule(name, role, feats, auth);

            BioDesign design = new BioDesign(name, auth);
            design.addPart(part);
            design.setModule(bMod);
            ObjectId id = persistor.save(design);
            
            return id;
        }
        else{
            System.out.println("Please check the parameter list passed to createPart, as it may be formatted incorrectly. Using createPart(name, author)...");
            return createPart(persistor, name, author);
        }
    }
}

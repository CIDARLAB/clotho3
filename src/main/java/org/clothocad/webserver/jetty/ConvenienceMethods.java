/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.webserver.jetty;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.List;
import org.ahocorasick.trie.Emit;
import org.ahocorasick.trie.Trie;
import org.ahocorasick.trie.Trie.TrieBuilder;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.model.Annotation;
import org.clothocad.model.BasicModule;
import org.clothocad.model.BioDesign;
import org.clothocad.model.Feature;
import org.clothocad.model.Module;
import org.clothocad.model.Parameter;
import org.clothocad.model.Part;
import org.clothocad.model.Person;
import org.clothocad.model.Sequence;

/**
 *
 * @author David
 */
public class ConvenienceMethods {

    //Base function
    public static ObjectId createPart(Persistor persistor, String name, String author) {
        Person auth = new Person(author);
        Part part = new Part(name, auth);

        BioDesign design = new BioDesign(name, auth);
        design.addPart(part);

        ObjectId id = persistor.save(design);

        return id;

    }

    //Optional list of parameters
    public static ObjectId createPart(Persistor persistor, String name, List<Parameter> parameters, String author) {
        Person auth = new Person(author);
        Part part = new Part(name, auth);

        BioDesign design = new BioDesign(name, auth);
        design.addPart(part);

        for (Parameter p : parameters) {
            design.addParameter(p);
        }

        ObjectId id = persistor.save(design);

        return id;
    }
    
    //Optional Map containing sequence and/or role
    public static ObjectId createPart(Persistor persistor, String name, Map<String, String> seqrole, String author) {
        if (seqrole.isEmpty()) {
            return createPart(persistor, name, author);
        }
        boolean bRole = false, bSeq = false;
        String role = "", sequence = "";

        //case insensitive check for sequence and/or role
        for (String field : seqrole.keySet()) {
            if (field.equalsIgnoreCase("role")) {
                bRole = true;
                role = seqrole.get(field);
            }

            if (field.equalsIgnoreCase("sequence")) {
                bSeq = true;
                sequence = seqrole.get(field);
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
        } else if (bSeq) {
            Sequence seq = new Sequence(name, sequence, auth);

            Part part = new Part(name, seq, auth);

            BioDesign design = new BioDesign(name, auth);
            design.addPart(part);

            ObjectId id = persistor.save(design);

            return id;

        } else if (bRole) {
            Part part = new Part(name, auth);
            Feature feat = new Feature(name, role, auth);
            
            Set<Feature> feats = new HashSet<>();
            feats.add(feat);

            Module bMod = new BasicModule(name, role, feats, auth);

            BioDesign design = new BioDesign(name, auth);
            design.addPart(part);
            design.setModule(bMod);
            ObjectId id = persistor.save(design);

            return id;
        } else {
            System.out.println("Please check the parameter list passed to createPart, as it may be formatted incorrectly. Using createPart(name, author)...");
            return createPart(persistor, name, author);
        }
    }

    //Both options
    public static ObjectId createPart(Persistor persistor, String name, Map<String, String> seqrole, List<Parameter> parameters, String author) {
        if (seqrole.isEmpty()) {
            return createPart(persistor, name, parameters, author);
        }
        boolean bRole = false, bSeq = false;
        String role = "", sequence = "";

        //case insensitive check for sequence and/or role
        for (String field : seqrole.keySet()) {
            if (field.equalsIgnoreCase("role")) {
                bRole = true;
                role = seqrole.get(field);
            }

            if (field.equalsIgnoreCase("sequence")) {
                bSeq = true;
                sequence = seqrole.get(field);
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

            for (Parameter p : parameters) {
                design.addParameter(p);
            }

            ObjectId id = persistor.save(design);

            return id;
        } else if (bSeq) {
            Sequence seq = new Sequence(name, sequence, auth);

            Part part = new Part(name, seq, auth);

            BioDesign design = new BioDesign(name, auth);
            design.addPart(part);

            for (Parameter p : parameters) {
                design.addParameter(p);
            }

            ObjectId id = persistor.save(design);

            return id;

        } else if (bRole) {
            Part part = new Part(name, auth);
            Feature feat = new Feature(name, role, auth);

            Set<Feature> feats = new HashSet<>();
            feats.add(feat);

            Module bMod = new BasicModule(name, role, feats, auth);

            BioDesign design = new BioDesign(name, auth);
            design.addPart(part);
            design.setModule(bMod);

            for (Parameter p : parameters) {
                design.addParameter(p);
            }

            ObjectId id = persistor.save(design);

            return id;
        } else {
            System.out.println("Please check the parameter list passed to createPart, as it may be formatted incorrectly. Using createPart(name, author)...");
            return createPart(persistor, name, author);
        }
    }
    
    //Variations of above that also have displayID option
    public static ObjectId createPart(Persistor persistor, String name, String displayID, String author){
        Person auth = new Person(author);
        Part part = new Part(name, auth);
        part.setDisplayID(displayID);

        BioDesign design = new BioDesign(name, auth);
        design.addPart(part);
        design.setDisplayID(displayID);

        ObjectId id = persistor.save(design);

        return id;
    }
    
    public static ObjectId createPart(Persistor persistor, String name, String displayID, List<Parameter> parameters, String author){
        Person auth = new Person(author);
        Part part = new Part(name, auth);
        part.setDisplayID(displayID);

        BioDesign design = new BioDesign(name, auth);
        design.addPart(part);
        design.setDisplayID(displayID);
        
        for(Parameter p : parameters)
        {
            design.addParameter(p);
        }

        ObjectId id = persistor.save(design);

        return id;
    }

    public static ObjectId createPart(Persistor persistor, String name, String displayID, Map<String, String> seqrole, String author){
        
        if (seqrole.isEmpty()) {
            return createPart(persistor, name, displayID, author);
        }
        boolean bRole = false, bSeq = false;
        String role = "", sequence = "";

        //case insensitive check for sequence and/or role
        for (String field : seqrole.keySet()) {
            if (field.equalsIgnoreCase("role")) {
                bRole = true;
                role = seqrole.get(field);
            }

            if (field.equalsIgnoreCase("sequence")) {
                bSeq = true;
                sequence = seqrole.get(field);
            }
        }

        Person auth = new Person(author);

        if (bSeq && bRole) {

            Sequence seq = new Sequence(name, sequence, auth);
            seq.setDisplayID(displayID);

            Part part = new Part(name, seq, auth);
            part.setDisplayID(displayID);
            
            Feature feat = new Feature(name, role, auth);
            feat.setSequence(seq);
            feat.setDisplayID(displayID);

            Annotation annotation = seq.createAnnotation(name, 1, sequence.length(), true, auth);
            annotation.setFeature(feat);

            Set<Feature> feats = new HashSet<>();
            feats.add(feat);

            Module bMod = new BasicModule(name, role, feats, auth);
            bMod.setDisplayID(displayID);
            
            BioDesign design = new BioDesign(name, auth);
            design.addPart(part);
            design.setModule(bMod);
            design.setDisplayID(displayID);
            ObjectId id = persistor.save(design);

            return id;
        } else if (bSeq) {
            Sequence seq = new Sequence(name, sequence, auth);

            Part part = new Part(name, seq, auth);
            part.setDisplayID(displayID);
            
            BioDesign design = new BioDesign(name, auth);
            design.addPart(part);
            design.setDisplayID(displayID);

            ObjectId id = persistor.save(design);

            return id;

        } else if (bRole) {
            Part part = new Part(name, auth);
            part.setDisplayID(displayID);
            
            Feature feat = new Feature(name, role, auth);
            feat.setDisplayID(displayID);

            Set<Feature> feats = new HashSet<>();
            feats.add(feat);

            Module bMod = new BasicModule(name, role, feats, auth);
            bMod.setDisplayID(displayID);
            
            BioDesign design = new BioDesign(name, auth);
            design.addPart(part);
            design.setModule(bMod);
            design.setDisplayID(displayID);
            
            ObjectId id = persistor.save(design);

            return id;
        } else {
            System.out.println("Please check the parameter list passed to createPart, as it may be formatted incorrectly. Using createPart(name, author)...");
            return createPart(persistor, name, author);
        }
    }
    
    public static ObjectId createPart(Persistor persistor, String name, String displayID, Map<String, String> seqrole, List<Parameter> parameters, String author){
        
        if (seqrole.isEmpty()) {
            return createPart(persistor, name, displayID, parameters, author);
        }
        boolean bRole = false, bSeq = false;
        String role = "", sequence = "";

        //case insensitive check for sequence and/or role
        for (String field : seqrole.keySet()) {
            if (field.equalsIgnoreCase("role")) {
                bRole = true;
                role = seqrole.get(field);
            }

            if (field.equalsIgnoreCase("sequence")) {
                bSeq = true;
                sequence = seqrole.get(field);
            }
        }

        Person auth = new Person(author);

        if (bSeq && bRole) {

            Sequence seq = new Sequence(name, sequence, auth);
            seq.setDisplayID(displayID);
            
            Part part = new Part(name, seq, auth);
            part.setDisplayID(displayID);
            
            Feature feat = new Feature(name, role, auth);
            feat.setSequence(seq);
            feat.setDisplayID(displayID);

            Annotation annotation = seq.createAnnotation(name, 1, sequence.length(), true, auth);
            annotation.setFeature(feat);

            Set<Feature> feats = new HashSet<>();
            feats.add(feat);

            Module bMod = new BasicModule(name, role, feats, auth);
            bMod.setDisplayID(displayID);

            BioDesign design = new BioDesign(name, auth);
            design.addPart(part);
            design.setModule(bMod);
            design.setDisplayID(displayID);

            for (Parameter p : parameters) {
                design.addParameter(p);
            }

            ObjectId id = persistor.save(design);

            return id;
        } else if (bSeq) {
            Sequence seq = new Sequence(name, sequence, auth);
            seq.setDisplayID(displayID);

            Part part = new Part(name, seq, auth);
            part.setDisplayID(displayID);
            
            BioDesign design = new BioDesign(name, auth);
            design.addPart(part);
            design.setDisplayID(displayID);

            for (Parameter p : parameters) {
                design.addParameter(p);
            }

            ObjectId id = persistor.save(design);

            return id;

        } else if (bRole) {
            Part part = new Part(name, auth);
            part.setDisplayID(displayID);
            
            Feature feat = new Feature(name, role, auth);
            feat.setDisplayID(displayID);

            Set<Feature> feats = new HashSet<>();
            feats.add(feat);

            Module bMod = new BasicModule(name, role, feats, auth);
            bMod.setDisplayID(displayID);

            BioDesign design = new BioDesign(name, auth);
            design.addPart(part);
            design.setModule(bMod);
            design.setDisplayID(displayID);
            
            for (Parameter p : parameters) {
                design.addParameter(p);
            }

            ObjectId id = persistor.save(design);

            return id;
        } else {
            System.out.println("Please check the parameter list passed to createPart, as it may be formatted incorrectly. Using createPart(name, author)...");
            return createPart(persistor, name, author);
        }
    }
    /*
        Create Device:
        
        No Role, No Sequence: 
            BioDesign
                ->subDesigns (createPart)
            Part(s)
                ->Assembly
                    ->subparts (Refers to Part objects in the subdesigns)
    
        With Role: (orange and blue)
            (BioDesign) -> BasicModule
                -> Feature
                    -> (if sequence) Sequence
                        -> Annotation
                            -> Feature (same one before Sequence)
    
        With Sequence: (red and blue)
            (Part) -> Sequence
                -> Annotations
                    -> To own feature
                    -> Features in subparts of the Assembly
        
     */
    
    //Base function
    public static ObjectId createDevice(Persistor persistor, String name, List<String> partIDs, String author, boolean createSeqFromParts) {

        Person auth = new Person(author);
        Part devPart = new Part(name, auth);
        devPart.createAssembly();

        BioDesign device = new BioDesign(name, auth);
        device.addPart(devPart);

        String seq = "";
        for (String id : partIDs) {
            ObjectId objId = new ObjectId(id);
            //partID refers to BioDesigns of parts
            BioDesign subDesign = persistor.get(BioDesign.class, objId);

            device.addSubDesign(subDesign);
            Set<Part> partSet = subDesign.getParts();

            for (Part p : partSet) {
                devPart.getAssemblies().get(0).addPart(p);
                if(createSeqFromParts)
                {
                    seq += p.getSequence().getSequence();
                }
            }
        }
        
        if(createSeqFromParts)
        {
            Sequence sequence = new Sequence(name, seq, auth);
            annotateMe(persistor, sequence, partIDs);
            devPart.setSequence(sequence);
        }

        ObjectId result = persistor.save(device);

        return result;
    }
    
    //Optional List of parameters
    public static ObjectId createDevice(Persistor persistor, String name, List<String> partIDs, List<Parameter> parameters, String author, boolean createSeqFromParts){
        
        Person auth = new Person(author);
        Part devPart = new Part(name, auth);
        devPart.createAssembly();

        BioDesign device = new BioDesign(name, auth);
        device.addPart(devPart);

        String seq = "";
        for (String id : partIDs) {
            ObjectId objId = new ObjectId(id);
            //partID refers to BioDesigns of parts
            BioDesign subDesign = persistor.get(BioDesign.class, objId);

            device.addSubDesign(subDesign);
            Set<Part> partSet = subDesign.getParts();

            for (Part p : partSet) {
                devPart.getAssemblies().get(0).addPart(p);
                if(createSeqFromParts)
                {
                    seq += p.getSequence().getSequence();
                }
            }
        }
        if(createSeqFromParts)
        {
            Sequence sequence = new Sequence(name, seq, auth);
            annotateMe(persistor, sequence, partIDs);
            devPart.setSequence(sequence);
        }
        
        for(Parameter p : parameters)
        {
            device.addParameter(p);
        }

        ObjectId result = persistor.save(device);

        return result;
    }

    //Optional Map containing sequence and/or role
    public static ObjectId createDevice(Persistor persistor, String name, List<String> partIDs, Map<String, String> seqrole, String author, boolean createSeqFromParts) {

        if (seqrole.isEmpty()) {
            return createDevice(persistor, name, partIDs, author, createSeqFromParts);
        } else {
            boolean bRole = false;
            boolean bSeq = createSeqFromParts;
            String role = "", sequence = "";
            //case insensitive check for parameters.
            for (String key : seqrole.keySet()) {
                if (key.equalsIgnoreCase("role")) {
                    bRole = true;
                    role = seqrole.get(key);
                }
                if (key.equalsIgnoreCase("sequence")) {
                    bSeq = true;
                    sequence = seqrole.get(key);
                }
            }

            Person auth = new Person(author);
            Part devPart = new Part(name, auth);
            devPart.createAssembly();

            BioDesign device = new BioDesign(name, auth);
            device.addPart(devPart);

            for (String id : partIDs) {
                ObjectId objId = new ObjectId(id);
                //partID refers to BioDesigns of parts
                BioDesign subDesign = persistor.get(BioDesign.class, objId);

                device.addSubDesign(subDesign);
                Set<Part> partSet = subDesign.getParts();

                for (Part p : partSet) {
                    devPart.getAssemblies().get(0).addPart(p);
                    //Even if we enable it, make sure you are NOT overwriting.
                    //Passing a sequence overrides the bool.
                    if(createSeqFromParts && sequence == "")
                    {
                        sequence += p.getSequence().getSequence();
                    }
                }
            }

            if (bSeq && bRole) {
                //Create objects
                Sequence seq = new Sequence(name, sequence, auth);
                Annotation seqAnno = seq.createAnnotation(name, 1, sequence.length(), true, auth);
                Feature annoFeat = new Feature(name, role, auth);

                //Chain objects together
                //Search and attach all other annotations from assembly to sequence
                annotateMe(persistor, seq, partIDs);
                //Attach seq to feature
                annoFeat.setSequence(seq);
                //Attach feature to annotation
                seqAnno.setFeature(annoFeat);
                //Attach annotation to sequence 
                seq.addAnnotation(seqAnno);
                //Attach sequence to part
                devPart.setSequence(seq);

                //Create a feature set with our feature
                HashSet<Feature> setFeat = new HashSet();
                setFeat.add(annoFeat);
                //Create a BasicModule with the feature set
                BasicModule bMod = new BasicModule(name, role, setFeat, auth);
                //Attach BasicModule to BioDesign
                device.setModule(bMod);
            } else if (bSeq) {

                //Create objects
                Sequence seq = new Sequence(name, sequence, auth);

                //Chain together objects
                //Annotate sequence with subfeatures from assembly and self
                annotateMe(persistor, seq, partIDs);
                
                //Attach sequence to part
                devPart.setSequence(seq);

            } else if (bRole) {

                //Create objects
                Feature annoFeat = new Feature(name, role, auth);
                HashSet<Feature> setFeat = new HashSet();
                setFeat.add(annoFeat);

                BasicModule bMod = new BasicModule(name, role, setFeat, auth);

                //Chain objects together
                //anno.setFeature(annoFeat);    //Commented out due to reasons above.
                device.setModule(bMod);
            }

            ObjectId result = persistor.save(device);

            return result;
        }
    }

    //Both options
    public static ObjectId createDevice(Persistor persistor, String name, List<String> partIDs, Map<String, String> seqrole, List<Parameter> parameters, String author, boolean createSeqFromParts){
        
        if (seqrole.isEmpty()) {
            return createDevice(persistor, name, partIDs, author, createSeqFromParts);
        } else {
            boolean bRole = false;
            boolean bSeq = createSeqFromParts;
            String role = "", sequence = "";
            //case insensitive check for parameters.
            for (String key : seqrole.keySet()) {
                if (key.equalsIgnoreCase("role")) {
                    bRole = true;
                    role = seqrole.get(key);
                }
                if (key.equalsIgnoreCase("sequence")) {
                    bSeq = true;
                    sequence = seqrole.get(key);
                }
            }

            Person auth = new Person(author);
            Part devPart = new Part(name, auth);
            devPart.createAssembly();

            BioDesign device = new BioDesign(name, auth);
            device.addPart(devPart);

            for (String id : partIDs) {
                ObjectId objId = new ObjectId(id);
                //partID refers to BioDesigns of parts
                BioDesign subDesign = persistor.get(BioDesign.class, objId);

                device.addSubDesign(subDesign);
                Set<Part> partSet = subDesign.getParts();

                for (Part p : partSet) {
                    devPart.getAssemblies().get(0).addPart(p);
                    //Even if we enable it, make sure you are NOT overwriting.
                    //Passing a sequence overrides the bool.                   
                    if(createSeqFromParts && sequence == "")
                    {
                        sequence += p.getSequence().getSequence();
                    }
                }
            }
            
            for(Parameter p : parameters){
                device.addParameter(p);
            }

            if (bSeq && bRole) {
                //Create objects
                Sequence seq = new Sequence(name, sequence, auth);
                Annotation seqAnno = seq.createAnnotation(name, 1, sequence.length(), true, auth);
                Feature annoFeat = new Feature(name, role, auth);

                //Chain objects together
                //Search and attach all other annotations from assembly to sequence
                annotateMe(persistor, seq, partIDs);
                //Attach seq to feature
                annoFeat.setSequence(seq);
                //Attach feature to annotation
                seqAnno.setFeature(annoFeat);
                //Attach annotation to sequence 
                seq.addAnnotation(seqAnno);
                //Attach sequence to part
                devPart.setSequence(seq);

                //Create a feature set with our feature
                HashSet<Feature> setFeat = new HashSet();
                setFeat.add(annoFeat);
                //Create a BasicModule with the feature set
                BasicModule bMod = new BasicModule(name, role, setFeat, auth);
                //Attach BasicModule to BioDesign
                device.setModule(bMod);
            } else if (bSeq) {

                //Create objects
                Sequence seq = new Sequence(name, sequence, auth);
                
                //Chain together objects
                //Annotate sequence with subfeatures from assembly and self
                annotateMe(persistor, seq, partIDs);
                
                //Attach sequence to part
                devPart.setSequence(seq);

            } else if (bRole) {

                //Create objects
                Feature annoFeat = new Feature(name, role, auth);
                HashSet<Feature> setFeat = new HashSet();
                setFeat.add(annoFeat);

                BasicModule bMod = new BasicModule(name, role, setFeat, auth);

                //Chain objects together
                //anno.setFeature(annoFeat);    //Commented out due to reasons above.
                device.setModule(bMod);
            }

            ObjectId result = persistor.save(device);

            return result;
        }
    }
    
    //Variations of above that also have the displayID option
    public static ObjectId createDevice(Persistor persistor, String name, List<String> partIDs, String displayID, String author, boolean createSeqFromParts){
        
        Person auth = new Person(author);
        Part devPart = new Part(name, auth);
        devPart.createAssembly();
        devPart.setDisplayID(displayID);

        BioDesign device = new BioDesign(name, auth);
        device.addPart(devPart);
        device.setDisplayID(displayID);

        String seq = "";
        for (String id : partIDs) {
            ObjectId objId = new ObjectId(id);
            //partID refers to BioDesigns of parts
            BioDesign subDesign = persistor.get(BioDesign.class, objId);

            device.addSubDesign(subDesign);
            Set<Part> partSet = subDesign.getParts();

            for (Part p : partSet) {
                devPart.getAssemblies().get(0).addPart(p);
                if(createSeqFromParts && p.getSequence() != null)
                {
                    seq += p.getSequence().getSequence();
                }
            }
        }
        
        if(createSeqFromParts)
        {
            Sequence sequence = new Sequence(name, seq, auth);
            annotateMe(persistor, sequence, partIDs);
            devPart.setSequence(sequence);
            sequence.setDisplayID(displayID);
        }

        ObjectId result = persistor.save(device);

        return result;
    }
    
    public static ObjectId createDevice(Persistor persistor, String name, List<String> partIDs, List<Parameter> parameters, String displayID, String author, boolean createSeqFromParts){
        
        Person auth = new Person(author);
        Part devPart = new Part(name, auth);
        devPart.createAssembly();
        devPart.setDisplayID(displayID);

        BioDesign device = new BioDesign(name, auth);
        device.addPart(devPart);
        device.setDisplayID(displayID);

        String seq = "";
        for (String id : partIDs) {
            ObjectId objId = new ObjectId(id);
            //partID refers to BioDesigns of parts
            BioDesign subDesign = persistor.get(BioDesign.class, objId);

            device.addSubDesign(subDesign);
            Set<Part> partSet = subDesign.getParts();

            for (Part p : partSet) {
                devPart.getAssemblies().get(0).addPart(p);
                if(createSeqFromParts)
                {
                    seq += p.getSequence().getSequence();
                }
            }
        }
        if(createSeqFromParts)
        {
            Sequence sequence = new Sequence(name, seq, auth);
            annotateMe(persistor, sequence, partIDs);
            devPart.setSequence(sequence);
            sequence.setDisplayID(displayID);
        }
        
        for(Parameter p : parameters)
        {
            device.addParameter(p);
        }

        ObjectId result = persistor.save(device);

        return result;
    }
    
    public static ObjectId createDevice(Persistor persistor, String name, List<String> partIDs, Map<String, String> seqrole, String displayID, String author, boolean createSeqFromParts){
        
        if (seqrole.isEmpty()) {
            return createDevice(persistor, name, partIDs, author, createSeqFromParts);
        } else {
            boolean bRole = false;
            boolean bSeq = createSeqFromParts;
            String role = "", sequence = "";
            //case insensitive check for parameters.
            for (String key : seqrole.keySet()) {
                if (key.equalsIgnoreCase("role")) {
                    bRole = true;
                    role = seqrole.get(key);
                }
                if (key.equalsIgnoreCase("sequence")) {
                    bSeq = true;
                    sequence = seqrole.get(key);
                }
            }

            Person auth = new Person(author);
            Part devPart = new Part(name, auth);
            devPart.createAssembly();
            devPart.setDisplayID(displayID);

            BioDesign device = new BioDesign(name, auth);
            device.addPart(devPart);
            device.setDisplayID(displayID);

            for (String id : partIDs) {
                ObjectId objId = new ObjectId(id);
                //partID refers to BioDesigns of parts
                BioDesign subDesign = persistor.get(BioDesign.class, objId);

                device.addSubDesign(subDesign);
                Set<Part> partSet = subDesign.getParts();

                for (Part p : partSet) {
                    devPart.getAssemblies().get(0).addPart(p);
                    //Even if we enable it, make sure you are NOT overwriting.
                    //Passing a sequence overrides the bool.
                    if(createSeqFromParts && sequence == "")
                    {
                        sequence += p.getSequence().getSequence();
                    }
                }
            }

            if (bSeq && bRole) {
                //Create objects
                Sequence seq = new Sequence(name, sequence, auth);
                Annotation seqAnno = seq.createAnnotation(name, 1, sequence.length(), true, auth);
                Feature annoFeat = new Feature(name, role, auth);

                //Chain objects together
                //Search and attach all other annotations from assembly to sequence
                annotateMe(persistor, seq, partIDs);
                //Attach seq to feature
                annoFeat.setSequence(seq);
                //Attach feature to annotation
                seqAnno.setFeature(annoFeat);
                //Attach annotation to sequence 
                seq.addAnnotation(seqAnno);
                //Attach sequence to part
                devPart.setSequence(seq);
                
                seq.setDisplayID(displayID);
                annoFeat.setDisplayID(displayID);

                //Create a feature set with our feature
                HashSet<Feature> setFeat = new HashSet();
                setFeat.add(annoFeat);
                //Create a BasicModule with the feature set
                BasicModule bMod = new BasicModule(name, role, setFeat, auth);
                bMod.setDisplayID(displayID);
                
                //Attach BasicModule to BioDesign
                device.setModule(bMod);
            } else if (bSeq) {

                //Create objects
                Sequence seq = new Sequence(name, sequence, auth);
                seq.setDisplayID(displayID);

                //Chain together objects
                //Annotate sequence with subfeatures from assembly and self
                annotateMe(persistor, seq, partIDs);
                
                //Attach sequence to part
                devPart.setSequence(seq);

            } else if (bRole) {

                //Create objects
                Feature annoFeat = new Feature(name, role, auth);
                HashSet<Feature> setFeat = new HashSet();
                setFeat.add(annoFeat);

                BasicModule bMod = new BasicModule(name, role, setFeat, auth);
                
                annoFeat.setDisplayID(displayID);
                bMod.setDisplayID(displayID);

                //Chain objects together
                //anno.setFeature(annoFeat);    //Commented out due to reasons above.
                device.setModule(bMod);
            }

            ObjectId result = persistor.save(device);

            return result;
        }
    }
    
    public static ObjectId createDevice(Persistor persistor, String name, List<String> partIDs, Map<String, String> seqrole, List<Parameter> parameters, String displayID, String author, boolean createSeqFromParts){
        
        if (seqrole.isEmpty()) {
            return createDevice(persistor, name, partIDs, author, createSeqFromParts);
        } else {
            boolean bRole = false;
            boolean bSeq = createSeqFromParts;
            String role = "", sequence = "";
            //case insensitive check for parameters.
            for (String key : seqrole.keySet()) {
                if (key.equalsIgnoreCase("role")) {
                    bRole = true;
                    role = seqrole.get(key);
                }
                if (key.equalsIgnoreCase("sequence")) {
                    bSeq = true;
                    sequence = seqrole.get(key);
                }
            }

            Person auth = new Person(author);
            Part devPart = new Part(name, auth);
            devPart.createAssembly();
            devPart.setDisplayID(displayID);

            BioDesign device = new BioDesign(name, auth);
            device.addPart(devPart);
            device.setDisplayID(displayID);

            for (String id : partIDs) {
                ObjectId objId = new ObjectId(id);
                //partID refers to BioDesigns of parts
                BioDesign subDesign = persistor.get(BioDesign.class, objId);

                device.addSubDesign(subDesign);
                Set<Part> partSet = subDesign.getParts();

                for (Part p : partSet) {
                    devPart.getAssemblies().get(0).addPart(p);
                    //Even if we enable it, make sure you are NOT overwriting.
                    //Passing a sequence overrides the bool.                   
                    if(createSeqFromParts && sequence == "")
                    {
                        sequence += p.getSequence().getSequence();
                    }
                }
            }
            
            for(Parameter p : parameters){
                device.addParameter(p);
            }

            if (bSeq && bRole) {
                //Create objects
                Sequence seq = new Sequence(name, sequence, auth);
                Annotation seqAnno = seq.createAnnotation(name, 1, sequence.length(), true, auth);
                Feature annoFeat = new Feature(name, role, auth);

                //Chain objects together
                //Search and attach all other annotations from assembly to sequence
                annotateMe(persistor, seq, partIDs);
                //Attach seq to feature
                annoFeat.setSequence(seq);
                //Attach feature to annotation
                seqAnno.setFeature(annoFeat);
                //Attach annotation to sequence 
                seq.addAnnotation(seqAnno);
                //Attach sequence to part
                devPart.setSequence(seq);

                //Create a feature set with our feature
                HashSet<Feature> setFeat = new HashSet();
                setFeat.add(annoFeat);
                //Create a BasicModule with the feature set
                BasicModule bMod = new BasicModule(name, role, setFeat, auth);
                //Attach BasicModule to BioDesign
                device.setModule(bMod);
                
                seq.setDisplayID(displayID);
                annoFeat.setDisplayID(displayID);
                bMod.setDisplayID(displayID);
            } else if (bSeq) {

                //Create objects
                Sequence seq = new Sequence(name, sequence, auth);
                seq.setDisplayID(displayID);
                
                //Chain together objects
                //Annotate sequence with subfeatures from assembly and self
                annotateMe(persistor, seq, partIDs);
                
                //Attach sequence to part
                devPart.setSequence(seq);

            } else if (bRole) {

                //Create objects
                Feature annoFeat = new Feature(name, role, auth);
                HashSet<Feature> setFeat = new HashSet();
                setFeat.add(annoFeat);

                BasicModule bMod = new BasicModule(name, role, setFeat, auth);

                annoFeat.setDisplayID(displayID);
                bMod.setDisplayID(displayID);
                
                //Chain objects together
                device.setModule(bMod);
            }

            ObjectId result = persistor.save(device);

            return result;
        }
    }

    //Scan the sequence for multiple string patterns, annotate it
    //Thank god someone invented the wheel (grep) already - Aho-Corasick Algorithm
    static void annotateMe(Persistor persistor, Sequence seq, List<String> partIDs) {
        HashMap<String, Part> partMap = new HashMap<>();

        TrieBuilder trieBuild = Trie.builder().removeOverlaps().caseInsensitive();
        //initialize list for string matching
        for (String s : partIDs) {
            ObjectId id = new ObjectId(s);
            BioDesign bd = persistor.get(BioDesign.class, id);
            Set<Part> parts = bd.getParts();
            for (Part p : parts) {
                if (p.getSequence() != null) {
                    if (p.getSequence().getSequence().isEmpty()) {
                        continue;
                    } else {
                        trieBuild.addKeyword(p.getSequence().getSequence());
                        partMap.put(p.getSequence().getSequence(), p);
                    }
                }
            }
        }
        Trie trie = trieBuild.build();

        Collection<Emit> results = trie.parseText(seq.getSequence());

        //Annotate seq based on results **DOUBLE CHECK THIS PROCESS WITH NIC
        //A lot of assumptions made, in practice I need to fully understand how 
        //data model will be used by Clotho/users.
        for (Emit em : results) {
            Part featPart = partMap.get(em.getKeyword());

            Feature feat = null;
            if (featPart.getSequence() == null) {
                continue;
            }
            else if(featPart.getSequence().getAnnotations() == null){
                continue;
            }

            System.out.println("AnnotateMe found : " + em.getKeyword());

            //Find part's annotation
            for (Annotation anno : featPart.getSequence().getAnnotations()) {
                if (anno.getFeature() != null) {
                    System.out.println("AnnotateMe annoFeat : " + anno.getFeature().getName());
                    if (anno.getFeature().getSequence() != null) {
                        //found feature that == part
                        System.out.println("Comparing with : " + anno.getFeature().getSequence().getSequence());

                        if (anno.getFeature().getSequence() == featPart.getSequence()) {
                            feat = anno.getFeature();
                            break;
                        }
                    }
                }
            }
            //if not null, annotate (should be most every time unless there's like no sequence or something).
            if (feat != null) {
                //add 1 to start and end because clotho starts counting from 1
                Annotation emAnno = seq.createAnnotation(featPart.getName(), em.getStart() + 1, em.getEnd() + 1, true, featPart.getAuthor());
                emAnno.setFeature(feat);
                seq.addAnnotation(emAnno);
            }
        }
    }
}

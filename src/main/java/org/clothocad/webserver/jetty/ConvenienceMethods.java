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
import java.util.Iterator;
import javax.persistence.EntityNotFoundException;
import org.ahocorasick.trie.Emit;
import org.ahocorasick.trie.Trie;
import org.ahocorasick.trie.Trie.TrieBuilder;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.model.Annotation;
import org.clothocad.model.Assembly;
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

    /*
    
    TO DO:
    Condense the create methods using the same format as the query methods,
    because as they are right now it's ugly and bad design.
    
    I started with only 2 iterations of each (the most basic forms), but both 
    times we needed to make revisions and add more parameters in a short time 
    span, the quick and dirty method was to just copy, paste, and adjust.
    
    __TLDR__ 
    Sorry, I know it's menial work but someone should condense the 
    create methods to match the query methods <3 - David T.
    
     */
 /*
    
        Create Methods
    
     */
    //////////////////////////
    //      Create Part     //
    //////////////////////////
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
    public static ObjectId createPart(Persistor persistor, String name, String displayID, String author) {
        Person auth = new Person(author);
        Part part = new Part(name, auth);
        part.setDisplayID(displayID);

        BioDesign design = new BioDesign(name, auth);
        design.addPart(part);
        design.setDisplayID(displayID);

        ObjectId id = persistor.save(design);

        return id;
    }

    public static ObjectId createPart(Persistor persistor, String name, String displayID, List<Parameter> parameters, String author) {
        Person auth = new Person(author);
        Part part = new Part(name, auth);
        part.setDisplayID(displayID);

        BioDesign design = new BioDesign(name, auth);
        design.addPart(part);
        design.setDisplayID(displayID);

        for (Parameter p : parameters) {
            design.addParameter(p);
        }

        ObjectId id = persistor.save(design);

        return id;
    }

    public static ObjectId createPart(Persistor persistor, String name, String displayID, Map<String, String> seqrole, String author) {

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

    public static ObjectId createPart(Persistor persistor, String name, String displayID, Map<String, String> seqrole, List<Parameter> parameters, String author) {

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

    //////////////////////////////
    //      Create Device       //
    //////////////////////////////
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
                if (createSeqFromParts) {
                    seq += p.getSequence().getSequence();
                }
            }
        }

        if (createSeqFromParts) {
            Sequence sequence = new Sequence(name, seq, auth);
            annotateMe(persistor, sequence, partIDs);
            devPart.setSequence(sequence);
        }

        ObjectId result = persistor.save(device);

        return result;
    }

    //Optional List of parameters
    public static ObjectId createDevice(Persistor persistor, String name, List<String> partIDs, List<Parameter> parameters, String author, boolean createSeqFromParts) {

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
                if (createSeqFromParts) {
                    seq += p.getSequence().getSequence();
                }
            }
        }
        if (createSeqFromParts) {
            Sequence sequence = new Sequence(name, seq, auth);
            annotateMe(persistor, sequence, partIDs);
            devPart.setSequence(sequence);
        }

        for (Parameter p : parameters) {
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
                    if (createSeqFromParts && sequence == "") {
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
    public static ObjectId createDevice(Persistor persistor, String name, List<String> partIDs, Map<String, String> seqrole, List<Parameter> parameters, String author, boolean createSeqFromParts) {

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
                    if (createSeqFromParts && sequence == "") {
                        sequence += p.getSequence().getSequence();
                    }
                }
            }

            for (Parameter p : parameters) {
                device.addParameter(p);
            }

            if (bSeq && bRole) {
                //Create objects
                Sequence seq = new Sequence(name, sequence, auth);
                Annotation seqAnno = seq.createAnnotation(name, 1, sequence.length(), true, auth);
                Feature annoFeat = new Feature(name, role, auth);
                
                //Create a feature set with our feature
                HashSet<Feature> setFeat = new HashSet();
                setFeat.add(annoFeat);
                //Create a BasicModule with the feature set
                BasicModule bMod = new BasicModule(name, role, setFeat, auth);
                //Attach BasicModule to BioDesign
                device.setModule(bMod);
                
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
    public static ObjectId createDevice(Persistor persistor, String name, String displayID, List<String> partIDs, String author, boolean createSeqFromParts) {

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
                if (createSeqFromParts && p.getSequence() != null) {
                    seq += p.getSequence().getSequence();
                }
            }
        }

        if (createSeqFromParts) {
            Sequence sequence = new Sequence(name, seq, auth);
            annotateMe(persistor, sequence, partIDs);
            devPart.setSequence(sequence);
            sequence.setDisplayID(displayID);
        }

        ObjectId result = persistor.save(device);

        return result;
    }

    public static ObjectId createDevice(Persistor persistor, String name, String displayID, List<String> partIDs, List<Parameter> parameters, String author, boolean createSeqFromParts) {

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
                if (createSeqFromParts) {
                    seq += p.getSequence().getSequence();
                }
            }
        }
        if (createSeqFromParts) {
            Sequence sequence = new Sequence(name, seq, auth);
            annotateMe(persistor, sequence, partIDs);
            devPart.setSequence(sequence);
            sequence.setDisplayID(displayID);
        }

        for (Parameter p : parameters) {
            device.addParameter(p);
        }

        ObjectId result = persistor.save(device);

        return result;
    }

    public static ObjectId createDevice(Persistor persistor, String name, String displayID, List<String> partIDs, Map<String, String> seqrole, String author, boolean createSeqFromParts) {

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
                    if (createSeqFromParts && sequence == "") {
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

    public static ObjectId createDevice(Persistor persistor, String name, String displayID, List<String> partIDs, Map<String, String> seqrole, List<Parameter> parameters, String author, boolean createSeqFromParts) {

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
                    if (createSeqFromParts && sequence == "") {
                        sequence += p.getSequence().getSequence();
                    }
                }
            }

            for (Parameter p : parameters) {
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

    /*
    
        Query Methods
    
     */
    //////////////////////
    //      Get Part    //
    //////////////////////
    //More often than not, exact = false will do the job.
    //I'd tiptoe around using exact = true unless you are absolutely certain you have exact fields (case sensitive)
    public static Map<String, Map<String, String>> getPart(Persistor persistor, Map<String, String> query) {
        return getPart(persistor, query.get("name"), query.get("displayID"), query.get("role"), query.get("sequence"), null, false);
    }

    public static Map<String, Map<String, String>> getPart(Persistor persistor, Map<String, String> query, List<Parameter> parameters, boolean exact) {
        return getPart(persistor, query.get("name"), query.get("displayID"), query.get("role"), query.get("sequence"), parameters, exact);
    }

    //Exact uses find(), broken for parameters
    public static Map<String, Map<String, String>> getPart(Persistor persistor, String name, String displayID, String role, String sequence, List<Parameter> parameters, boolean exact) {
        HashMap<String, Object> query = new HashMap<>();
        Iterable<ObjBase> nameDisplayParameterList = null;
        Iterable<ObjBase> sequenceList = null;
        Iterable<ObjBase> roleList = null;

        //nameDisplayParameter list
        if (displayID != null) {
            query.put("displayID", displayID);
        }
        if (name != null) {
            query.put("name", name);
        }
        if (parameters != null) {
            if (parameters.size() > 1) {
                String valuePrefix = "{$regex: '";
                String valueSuffix = "', $options: 'i'}";
                for (Parameter p : parameters) {
                    query.put("parameters multiMatch " + p.getName(),
                            "{name: " + valuePrefix + p.getName() + valueSuffix
                            + ", variable: " + valuePrefix + p.getVariable() + valueSuffix
                            + ", units: " + valuePrefix + p.getUnits() + valueSuffix + "}"
                    );
                }
            } else {
                query.put("parameters.name", parameters.get(0).getName());
                query.put("parameters.variable", parameters.get(0).getVariable());
                query.put("parameters.units", parameters.get(0).getUnits());
            }
        }
        if (!query.isEmpty()) {
            nameDisplayParameterList = (exact) ? persistor.find(query) : persistor.findRegex(query);
            query = new HashMap<>();
        }

        //sequenceList
        if (sequence != null) {
            if (name != null) {
                query.put("name", name);
            }
            query.put("sequence", sequence);
            sequenceList = (exact) ? persistor.find(query) : persistor.findRegex(query);
            query = new HashMap<>();
        }

        //roleList
        if (role != null) {
            if (name != null) {
                query.put("name", name);
            }
            query.put("role", role);
            roleList = (exact) ? persistor.find(query) : persistor.findRegex(query);
            query = new HashMap<>();
        }

        //Retrieve objects that matched
        Map<String, Map<String, String>> returnList = new HashMap<>();

        /*
        Assumptions:
        
        nameDisplayParameterList SHOULD have found it. sequenceList and roleList
        therefore only have to append role and sequence fields.
        
        If for some reason it wasn't found by nameDisplayParameterList, assume
        that parameters did not match/don't exist. sequenceList will create an
        entry with its own name, displayID (if exists), and sequence. roleList
        will then append the role field.
        
        If for some reason it wasn't found by any of the above, roleList will
        create an entry with its own name, displayID, and role.
         */
        if (nameDisplayParameterList != null) {
            for (ObjBase each : nameDisplayParameterList) {
                BioDesign bd = persistor.get(BioDesign.class, each.getId());
                if (bd != null) {
                    if (!returnList.containsKey(bd.getName())) {
                        Map<String, String> insert = new HashMap<>();
                        insert.put("name", bd.getName());
                        insert.put("id", bd.getId().toString());
                        if (bd.getDisplayID() != null) {
                            insert.put("displayID", bd.getDisplayID());
                        }
                        if (bd.getParameters() != null) {
                            String build = "[";
                            for (Parameter p : bd.getParameters()) {
                                build += "{name:'" + p.getName()
                                        + "', value:" + p.getValue()
                                        + ", variable:'" + p.getVariable()
                                        + "', units:'" + p.getUnits() + "'},";
                            }
                            insert.put("parameters", build + "]");
                        }
                        returnList.put(bd.getName(), insert);
                    }
                }
            }
        }
        if (sequenceList != null) {
            for (ObjBase each : sequenceList) {
                Sequence s = persistor.get(Sequence.class, each.getId());
                if (s != null) {
                    if (returnList.containsKey(s.getName())) {
                        returnList.get(s.getName()).put("sequence", s.getSequence());
                    } else {
                        HashMap<String, String> insert = new HashMap<>();
                        insert.put("name", s.getName());
                        if (s.getDisplayID() != null) {
                            insert.put("displayID", s.getDisplayID());
                        }
                        insert.put("sequence", s.getSequence());
                        returnList.put(s.getName(), insert);
                    }
                }
            }
        }
        if (roleList != null) {
            for (ObjBase each : roleList) {
                if (persistor.get(each.getId()).getClass().equals(Feature.class)) {
                    Feature f = (Feature) persistor.get(each.getId());
                    if (f != null) {
                        if (returnList.containsKey(f.getName())) {
                            returnList.get(f.getName()).put("role", f.getRole());
                        } else {
                            HashMap<String, String> insert = new HashMap<>();
                            insert.put("name", f.getName());
                            if (f.getDisplayID() != null) {
                                insert.put("displayID", f.getDisplayID());
                            }
                            insert.put("role", f.getRole());
                            returnList.put(f.getName(), insert);
                        }
                    }
                } else {
                    Module m = (Module) persistor.get(each.getId());
                    if (m != null) {
                        if (returnList.containsKey(m.getName())) {
                            returnList.get(m.getName()).put("role", m.getRole());
                        } else {
                            HashMap<String, String> insert = new HashMap<>();
                            insert.put("name", m.getName());
                            if (m.getDisplayID() != null) {
                                insert.put("displayID", m.getDisplayID());
                            }
                            insert.put("role", m.getRole());
                            returnList.put(m.getName(), insert);
                        }
                    }
                }
            }
        }
        return returnList;
    }

    public static List<String> getPartID(Persistor persistor, Map<String, String> query, List<Parameter> parameters, boolean exact) {
        Map<String, Map<String, String>> res = getPart(persistor, query, parameters, exact);
        ArrayList<String> ids = new ArrayList<>();

        for (String key : res.keySet()) {
            if (res.get(key).containsKey("id")) {
                ids.add(res.get(key).get("id"));
            }
        }

        return ids;
    }
    //////////////////////////
    //      Get Device      //
    //////////////////////////
    //More often than not, exact = false will do the job.
    //I'd tiptoe around using exact = true unless you are absolutely certain you have exact fields (case sensitive)

    public static Map<String, Map<String, String>> getDevice(Persistor persistor, Map<String, String> query) {
        return getDevice(persistor, query.get("name"), query.get("displayID"), query.get("role"), query.get("sequence"), null, null, false);
    }

    public static Map<String, Map<String, String>> getDevice(Persistor persistor, Map<String, String> query, Map<String, List> subObjects, boolean exact) {
        return getDevice(persistor, query.get("name"), query.get("displayID"), query.get("role"), query.get("sequence"), (List<String>) subObjects.get("parts"), (List<Parameter>) subObjects.get("parameters"), exact);
    }

    public static Map<String, Map<String, String>> getDevice(Persistor persistor, List<String> parts) {
        return getDevice(persistor, null, null, null, null, parts, null, false);
    }

    //[String displayID], [String name], [String role], [String sequence], [Part[] parts]
    public static Map<String, Map<String, String>> getDevice(Persistor persistor, String name, String displayID, String role, String sequence, List<String> parts, List<Parameter> parameters, boolean exact) {

        HashMap<String, Object> query = new HashMap<>();
        Iterable<ObjBase> nameDisplayParameterList = null;
        Iterable<ObjBase> sequenceList = null;
        Iterable<ObjBase> roleList = null;

        //nameDisplayParameter list
        if (displayID != null) {
            query.put("displayID", displayID);
        }
        if (name != null) {
            query.put("name", name);
        }
        if (parameters != null) {
            if (parameters.size() > 1) {
                String valuePrefix = "{$regex: '";
                String valueSuffix = "', $options: 'i'}";
                for (Parameter p : parameters) {
                    query.put("parameters multiMatch " + p.getName(),
                            "{name: " + valuePrefix + p.getName() + valueSuffix
                            + ", variable: " + valuePrefix + p.getVariable() + valueSuffix
                            + ", units: " + valuePrefix + p.getUnits() + valueSuffix + "}"
                    );
                }
            } else {
                query.put("parameters.name", parameters.get(0).getName());
                query.put("parameters.variable", parameters.get(0).getVariable());
                query.put("parameters.units", parameters.get(0).getUnits());
            }
        }

        if (parts != null) {

            for (String id : parts) {
                query.put("parts multiMatch", id);
            }

        }

        if (!query.isEmpty()) {
            nameDisplayParameterList = (exact) ? persistor.find(query) : persistor.findRegex(query);
            query = new HashMap<>();
        }

        //sequenceList
        if (sequence != null) {
            if (name != null) {
                query.put("name", name);
            }
            query.put("sequence", sequence);
            sequenceList = (exact) ? persistor.find(query) : persistor.findRegex(query);
            query = new HashMap<>();
        }

        //roleList
        if (role != null) {
            if (name != null) {
                query.put("name", name);
            }
            query.put("role", role);
            roleList = (exact) ? persistor.find(query) : persistor.findRegex(query);
            query = new HashMap<>();
        }

        //Retrieve objects that matched
        Map<String, Map<String, String>> returnList = new HashMap<>();

        /*
        Assumptions:
        
        nameDisplayParameterList SHOULD have found it. sequenceList and roleList
        therefore only have to append role and sequence fields.
        
        If for some reason it wasn't found by nameDisplayParameterList, assume
        that parameters did not match/don't exist. sequenceList will create an
        entry with its own name, displayID (if exists), and sequence. roleList
        will then append the role field.
        
        If for some reason it wasn't found by any of the above, roleList will
        create an entry with its own name, displayID, and role.
         */
        if (nameDisplayParameterList != null) {
            for (ObjBase each : nameDisplayParameterList) {

                BioDesign bd = persistor.get(BioDesign.class, each.getId());
                if (bd != null) {
                    if (!returnList.containsKey(bd.getName())) {
                        Map<String, String> insert = new HashMap<>();
                        insert.put("name", bd.getName());
                        insert.put("id", bd.getId().toString());
                        if (bd.getDisplayID() != null) {
                            insert.put("displayID", bd.getDisplayID());
                        }
                        if (bd.getParameters() != null) {
                            String build = "[";
                            for (Parameter p : bd.getParameters()) {
                                build += "{name:'" + p.getName()
                                        + "', value:" + p.getValue()
                                        + ", variable:'" + p.getVariable()
                                        + "', units:'" + p.getUnits() + "'},";
                            }
                            insert.put("parameters", build + "]");
                        }
                        //Get Device also needs to list parts.
                        if (bd.getParts() != null) {
                            System.out.println("Query found parts for " + bd.getName());
                            String pbuild = "[";
                            for (Part p : bd.getParts()) {
                                pbuild += "{id: '" + p.getId() + "', name:'" + p.getName()
                                        + "', displayID:'" + p.getDisplayID();
                                if (p.getSequence() != null) {
                                    if (p.getSequence().getAnnotations() != null) {
                                        for (Annotation anno : p.getSequence().getAnnotations()) {
                                            if (anno.getFeature().getName().equalsIgnoreCase(p.getName())) {
                                                pbuild += "', role:'" + anno.getFeature().getName();
                                                break;
                                            }
                                        }
                                    }
                                    pbuild += "', sequence:'" + p.getSequence().getSequence();
                                }
                                pbuild += "'},";
                            }
                            insert.put("parts", pbuild + "]");
                        }
                        System.out.println(insert.values());
                        returnList.put(bd.getName(), insert);
                    }
                }
            }
        }
        if (sequenceList != null) {
            for (ObjBase each : sequenceList) {
                Sequence s = persistor.get(Sequence.class, each.getId());
                if (s != null) {
                    if (returnList.containsKey(s.getName())) {
                        returnList.get(s.getName()).put("sequence", s.getSequence());
                    } else {
                        HashMap<String, String> insert = new HashMap<>();
                        insert.put("name", s.getName());
                        if (s.getDisplayID() != null) {
                            insert.put("displayID", s.getDisplayID());
                        }
                        insert.put("sequence", s.getSequence());
                        returnList.put(s.getName(), insert);
                    }
                }
            }
        }
        if (roleList != null) {
            for (ObjBase each : roleList) {
                if (persistor.get(each.getId()).getClass().equals(Feature.class)) {
                    Feature f = (Feature) persistor.get(each.getId());
                    if (f != null) {
                        if (returnList.containsKey(f.getName())) {
                            returnList.get(f.getName()).put("role", f.getRole());
                        } else {
                            HashMap<String, String> insert = new HashMap<>();
                            insert.put("name", f.getName());
                            if (f.getDisplayID() != null) {
                                insert.put("displayID", f.getDisplayID());
                            }
                            insert.put("role", f.getRole());
                            returnList.put(f.getName(), insert);
                        }
                    }
                } else {
                    Module m = (Module) persistor.get(each.getId());
                    if (m != null) {
                        if (returnList.containsKey(m.getName())) {
                            returnList.get(m.getName()).put("role", m.getRole());
                        } else {
                            HashMap<String, String> insert = new HashMap<>();
                            insert.put("name", m.getName());
                            if (m.getDisplayID() != null) {
                                insert.put("displayID", m.getDisplayID());
                            }
                            insert.put("role", m.getRole());
                            returnList.put(m.getName(), insert);
                        }
                    }
                }
            }
        }
        return returnList;
    }

    public static List<String> getDeviceID(Persistor persistor, Map<String, String> query, Map<String, List> subObjects, boolean exact) {
        Map<String, Map<String, String>> res = getDevice(persistor, query, subObjects, exact);
        ArrayList<String> ids = new ArrayList<>();

        for (String key : res.keySet()) {
            if (res.get(key).containsKey("id")) {
                ids.add(res.get(key).get("id"));
            }
        }

        return ids;

    }

    //////////////////////
    //      Delete      //
    //////////////////////
    public static boolean delete(Persistor persistor, ObjectId id, boolean bDevice) {
        try {
            BioDesign top = persistor.get(BioDesign.class, id);

            deleteReferences(persistor, id);

            if (top.getSubDesigns() != null && bDevice) {
                for (BioDesign bd : top.getSubDesigns()) {
                    boolean bDeleted = delete(persistor, bd.getId(), bDevice);
                    System.out.println(bd.getId() + " : " + bDeleted);
                }
            }

            ArrayList<ObjBase> deleteList = new ArrayList<>();
            deleteList.add(top);
            if (top.getModule() != null) {
                deleteList.add(top.getModule());
            }

            for (Part p : top.getParts()) {
                if (p.getId() != null) {
                    deleteList.add(p);
                }
                for (Annotation anno : p.getSequence().getAnnotations()) {
                    if (anno.getId() != null) {
                        deleteList.add(anno);

                        if (anno.getFeature() != null) {
                            deleteList.add(anno.getFeature());
                        }
                    }
                }
                if (p.getSequence() != null) {
                    deleteList.add(p.getSequence());
                }
            }

            persistor.delete(deleteList);

            return true;
        } catch (Exception e) {

            System.out.println("[ERROR] Convenience Delete : " + e.toString());
            e.printStackTrace();
            return false;
        }
    }

    //////////////////////////////////
    //      Utility Functions       //
    //////////////////////////////////
    //Scan a sequence for multiple string patterns, annotate it
    //Thank god someone invented the wheel (grep) already - Aho-Corasick Algorithm
    static void annotateMe(Persistor persistor, Sequence seq, List<String> partIDs) {
        HashMap<String, Part> partMap = new HashMap<>();

        TrieBuilder trieBuild = Trie.builder().removeOverlaps().caseInsensitive();
        //initialize list for string matching
        for (String s : partIDs) {
            ObjectId id = new ObjectId(s);
            BioDesign bd = persistor.get(BioDesign.class, id);
            if (bd != null) {
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
            } else if (featPart.getSequence().getAnnotations() == null) {
                continue;
            }
//            System.out.println("AnnotateMe found : " + em.getKeyword());

            //Find part's annotation
            for (Annotation anno : featPart.getSequence().getAnnotations()) {
                if (anno.getFeature() != null) {
//                    System.out.println("AnnotateMe annoFeat : " + anno.getFeature().getName());
                    if (anno.getFeature().getSequence() != null) {
                        //found feature that == part
//                        System.out.println("Comparing with : " + anno.getFeature().getSequence().getSequence());

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

    public static void deleteReferences(Persistor persistor, ObjectId bioId) {
        ObjectId partId = null;

        //Get the a list of devices that include partId in their parts list
        ArrayList<String> query = new ArrayList<>();
        query.add(bioId.getValue());
        Map<String, Map<String, String>> ref = getDevice(persistor, query);
        for (String key : ref.keySet()) {
            System.out.println("Deleting reference from " + key);

            String id = ref.get(key).get("id");
            BioDesign obj = persistor.get(BioDesign.class, new ObjectId(id));

            //Remove from the subDesigns
            //list of parts, the annotations, etc.
            for (BioDesign each : obj.getSubDesigns()) {
                if (each.getId().getValue().equalsIgnoreCase(bioId.getValue())) {
                    for (Part p : obj.getParts()) {
                        if (p.getName().equalsIgnoreCase(each.getName())) {

                            System.out.println("Deleting subdesign " + bioId.getValue());
                            partId = p.getId();
                            break;
                        }
                    }
                    obj.getSubDesigns().remove(each);
                    break;
                }
            }
            //If part exists, remove from Assemblies
            if (partId != null) {
                for (Part p : obj.getParts()) {
                    List<Assembly> newA = new ArrayList<>();
                    for (Assembly a : p.getAssemblies()) {
                        for (Part aPart : a.getParts()) {
                            if (aPart.getId().getValue().equalsIgnoreCase(partId.getValue())) {
                                System.out.println("Removing from assembly part " + partId.getValue());
                                a.getParts().remove(aPart);
                                break;
                            }
                        }
                    }
                }
            }

            //Overwrite existing Device with smaller Assemblies/Sub-Designs
            System.out.println("Overwriting " + key);
            persistor.save(obj, true);
        }

    }

    /**
     * ********************
     *
     * Author: Jerome
     *
     * *******************
     */
    private static void updateBioDesign(
            Persistor persistor,
            ObjectId obj,
            String displayID,
            String name,
            List<Parameter> parameters,
            List<String> subPartIds,
            Map<String, String> seqrole,
            boolean createSeqFromParts) {

        BioDesign bio = persistor.get(BioDesign.class, obj);

        // Update biodesign's displayID
        if (displayID != null) {
            bio.setDisplayID(displayID);
        }

        // Update biodesign's name
        if (name != null) {
            bio.setName(name);
        }

        // Update biodesign's parameters
        if (parameters != null) {
            // There exists a Set<Parameters>
            if (bio.getParameters() != null) {
                bio.getParameters().clear();
                for (Parameter parameter : parameters) {
                    bio.getParameters().add(parameter);
                }
            } else { //No Set<Parameters>; .addParameter creates new hashset
                for (Parameter parameter : parameters) {
                    bio.addParameter(parameter);
                }
            }
        }

        // Update biodesign's subPartIds
        if (subPartIds != null) {
            // Check for null pointer. Should always be true.
            
            if (bio.getParts().size() > 0) {

                String sequenceString = "";

                Part[] partArray = bio.getParts().toArray(new Part[bio.getParts().size()]);
                Part part = partArray[0];   // Nominal part

                // Clear if there is an existing Set<BioDesign>
                if (bio.getSubDesigns() != null) {
                    bio.getSubDesigns().clear();
                }

                // Clear if there is an existing assembly
                if (part.getAssemblies() != null) {
                    part.getAssemblies().clear();
                }

                Assembly assembly = part.createAssembly();

                // Add to assembly and subdesigns
                for (String subPartId : subPartIds) {
                    ObjectId id = new ObjectId(subPartId);

                    BioDesign subPart = persistor.get(BioDesign.class, id);
                    Part[] subPartArray = subPart.getParts().toArray(new Part[subPart.getParts().size()]);

                    assembly.addPart(subPartArray[0]);
                    bio.addSubDesign(subPart);

                    sequenceString += subPartArray[0].getSequence().getSequence();
                }

                if (createSeqFromParts) {

                    Sequence seq;
                    Annotation seqAnno;
                    Set<Feature> featSet = new HashSet<>();

                    if (part.getSequence() == null) {
                        seq = new Sequence(name, sequenceString, bio.getAuthor());
                        seqAnno = seq.createAnnotation(name, 1, sequenceString.length(), true, bio.getAuthor());
                        part.setSequence(seq);

                    } else {
                        seq = part.getSequence();
                        part.getSequence().setSequence(sequenceString);
                        Annotation[] annoArray = part.getSequence().getAnnotations().toArray(new Annotation[part.getSequence().getAnnotations().size()]);
                        seqAnno = annoArray[0];
                        seqAnno.setEnd(sequenceString.length());
                    }

                    seq.getAnnotations().clear();
                    annotateMe(persistor, seq, subPartIds);
                    
                    
                    if (bio.getModule() != null) {
                        if (seqAnno.getFeature() != null) {
                            featSet.add(seqAnno.getFeature());
                        }
                    }
                }
            } else {
                System.out.println("There are no parts in this BioDevice.");
            }
        }

        // Update biodesign's seqrole  
        if (seqrole != null) {
            if (bio.getParts().size() > 0) {
                Part[] partArray = bio.getParts().toArray(new Part[bio.getParts().size()]);
                Part part = partArray[0];

                boolean bRole = false, bSeq = false;
                String roleString = "", sequenceString = "";

                //case insensitive check for sequence and/or role
                for (String field : seqrole.keySet()) {
                    if (field.equalsIgnoreCase("role")) {
                        bRole = true;
                        roleString = seqrole.get(field);
                    }

                    if (field.equalsIgnoreCase("sequence")) {
                        bSeq = true;
                        sequenceString = seqrole.get(field);
                    }
                }

                if (bSeq && bRole) {
                    Sequence sequence = part.getSequence();
                    sequence.setSequence(sequenceString);

                    Annotation[] annos = sequence.getAnnotations().toArray(new Annotation[sequence.getAnnotations().size()]);
                    Feature feature = annos[0].getFeature();

                    annos[0].setEnd(sequenceString.length());
                    annos[0].getFeature().setRole(roleString);
                    
                    if (bio.getModule() != null) {
                        bio.getModule().setRole(roleString);
                    } else {
                        Set<Feature> setFeat = new HashSet<>();
                        setFeat.add(feature);
                        BasicModule module = new BasicModule(name, roleString, setFeat, bio.getAuthor());
                        bio.setModule(module);
                        bio.getModule().setRole(roleString);
                    }
                    
                } else if (bSeq){

                    Sequence sequence = part.getSequence();
                    sequence.setSequence(sequenceString);

                    Annotation[] annos = sequence.getAnnotations().toArray(new Annotation[sequence.getAnnotations().size()]);
                    annos[0].setEnd(sequenceString.length());

                } else if (bRole) {
                    Sequence sequence = part.getSequence();
                    if (sequence == null) {
                        System.out.println("Old Sequence is NULL. Role will not be updated.");
                        return;
                    }

                    Annotation[] annos = sequence.getAnnotations().toArray(new Annotation[sequence.getAnnotations().size()]);
                    Feature feature = annos[0].getFeature();
                    feature.setRole(roleString);
                    if (bio.getModule() != null) {
                        bio.getModule().setRole(roleString);
                    } else {
                        Set<Feature> featSet = new HashSet<>();
                        featSet.add(feature);
                        BasicModule module = new BasicModule(name, roleString, featSet, bio.getAuthor());
                        bio.setModule(module);
                    }

                } else {
                    System.out.println("seqRole formatted incorrectly.");
                }

            }
        }

        persistor.save(bio, true);
    }

    public static void updatePart(
            Persistor persistor,
            ObjectId obj,
            String displayID,
            String name,
            List<Parameter> parameters,
            Map<String, String> seqrole) {

        updateBioDesign(persistor, obj, displayID, name, parameters, null, seqrole, false);
    }

    public static void updateDevice(
            Persistor persistor,
            ObjectId obj,
            String displayID,
            String name,
            List<Parameter> parameters,
            List<String> subPartIds,
            Map<String, String> seqrole,
            boolean createSeqFromParts) {

        updateBioDesign(persistor, obj, displayID, name, parameters, subPartIds, seqrole, createSeqFromParts);
    }
}

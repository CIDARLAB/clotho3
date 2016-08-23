package org.clothocad.model;

import java.util.HashMap;
import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.Reference;
import org.clothocad.core.persistence.annotations.ReferenceCollection;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import org.json.simple.JSONObject;

/**
*
* @author Nicholas Roehner
* @author mardian
*/
@NoArgsConstructor
public class BioDesign extends SharableObjBase {

    @Getter
    @Setter
    @Reference
    protected Module module;

    @Getter
    protected Set<Parameter> parameters;

    @Getter
    @Setter
    @ReferenceCollection
    protected Set<Part> parts;

    @Getter
    @Setter
    @ReferenceCollection
    protected Set<Polynucleotide> polynucleotides;

    @Getter
    @Setter
    @ReferenceCollection
    protected Set<Strain> strains;

    @Getter
    @Setter
    @ReferenceCollection
    protected Set<Medium> media;

    @Getter
    @Setter
    @ReferenceCollection
    protected Set<BioDesign> subDesigns;

    @Getter
    @Setter
    @Reference
    protected BioDesign parentDesign;

    public BioDesign(String name, String description, Person author) {
        super(name, author, description);
    }

    public BioDesign(String name, String description, BasicModule module, Person author) {
        super(name, author, description);
        this.module = module;
    }

    public BioDesign(String name, String description, CompositeModule module, Person author) {
        super(name, author, description);
        this.module = module;
    }

    public BioDesign(String name, Person author) {
        super(name, author);
    }

    public Parameter createParameter(double value, Variable variable) {
        Parameter parameter = new Parameter(value, variable);
        addParameter(parameter);
        return parameter;
    }
    
    public void addParameter(Parameter parameter) {
    	if (parameters == null) {
    		parameters = new HashSet<Parameter>();
    	}
    	parameters.add(parameter);
    }
    
    public void addPart(Part part) {
    	if (parts == null) {
    		parts = new HashSet<Part>();
    	}
    	parts.add(part);
    }
    
    public void addPolynucleotide(Polynucleotide polynucleotide) {
    	if (polynucleotides == null) {
    		polynucleotides = new HashSet<Polynucleotide>();
    	}
    	polynucleotides.add(polynucleotide);
    }
    
    public void addStrain(Strain strain) {
    	if (strains == null) {
    		strains = new HashSet<Strain>();
    	}
    	strains.add(strain);
    }
    
    public void addMedium(Medium medium) {
    	if (media == null) {
    		media = new HashSet<Medium>();
    	}
    	media.add(medium);
    }
    
    public void addSubDesign(BioDesign subDesign) {
    	if (subDesigns == null) {
    		subDesigns = new HashSet<BioDesign>();
    	}
    	subDesigns.add(subDesign);
    }

    public Map getMap(){
        Map map = new HashMap();
        map.put("name", this.getName());
        map.put("author", this.getAuthor().getName());
        map.put("description", this.getDescription());
        if(this.module!=null)
            map.put("module", this.module.getName());
        if(this.parameters!=null)
            map.put("parameters", this.parameters.size());
        if(this.parts!=null)
            map.put("parts", this.parts.size());
        if(this.polynucleotides!=null)
            map.put("polynucleotides", this.polynucleotides.size());
        if(this.strains!=null)
            map.put("strains", this.strains.size());
        if(this.media!=null)
            map.put("media", this.media.size());
        if(this.subDesigns!=null)
            map.put("subDesigns", this.subDesigns.size());
        return map;
    }
    
    public JSONObject getJSON(){
        JSONObject obj = new JSONObject();
        obj.put("name", this.getName());
        obj.put("author", this.getAuthor().getName());
        obj.put("description", this.getDescription());
        if(this.module!=null)
            obj.put("module", this.module.getName());
        if(this.parameters!=null)
            obj.put("parameters", this.parameters.size());
        if(this.parts!=null)
            obj.put("parts", this.parts.size());
        if(this.polynucleotides!=null)
            obj.put("polynucleotides", this.polynucleotides.size());
        if(this.strains!=null)
            obj.put("strains", this.strains.size());
        if(this.media!=null)
            obj.put("media", this.media.size());
        if(this.subDesigns!=null)
            obj.put("subDesigns", this.subDesigns.size());
        return obj;
    }
    
    /*public String toString(){
        String str = "";
        str += "Name : " + this.getName() + "\n";
        str += "Author : " + this.getAuthor().getName() + "\n";
        str += "Description : " + this.getDescription() + "\n";
        str += "Sequence : " + this.sequence.getSequence() + "\n";
        str += "GenBank ID : " + this.genbankId + "\n";
        str += "SwissProt ID : " + this.swissProtId + "\n";
        str += "PDB ID : " + this.PDBId + "\n";
        str += "Risk Group : " + this.riskGroup + "\n";
        str += "Feature Role : " + this.role.toString();
        if (this.parentFeature!=null)
            str += "\nParent Feature : " + this.parentFeature.getName();
        return str;
    }*/
}

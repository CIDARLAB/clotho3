package org.clothocad.model;

import java.util.HashMap;
import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.ReferenceCollection;
import org.clothocad.core.persistence.annotations.Reference;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import javax.validation.constraints.NotNull;
import javax.validation.constraints.Size;
import org.json.simple.JSONObject;

/**
*
* @author Nicholas Roehner
* @author mardian
*/
@NoArgsConstructor
public class ExperimentalDesign extends SharableObjBase {

    @NotNull
    @Size(min=1)
    @Getter
    @Setter
    @ReferenceCollection
    protected Set<Variable> responseVariables;

    @NotNull
    @Size(min=1)
    @Getter
    @Setter
    @ReferenceCollection
    protected Set<Variable> controlledVariables;

    @Getter
    protected Set<ExperimentalCondition> experimentalConditions;

    @Getter
    @Setter
    @Reference
    protected BioDesign bioDesign;

    @Getter
    @Setter
    @ReferenceCollection
    protected Set<ExperimentalDesign> subDesigns;
    
    @Getter
    @Setter
    @Reference
    protected ExperimentalDesign parentDesign;
    
    @Getter
    @Setter
    @Reference
    protected int site;
    
    @Getter
    @Setter
    @Reference
    protected String notes;
    
    public ExperimentalDesign(String name, Set<Variable> responseVariables, Set<Variable> controlledVariables,
            Person author) {
        super(name, author);
        this.responseVariables = responseVariables;
        this.controlledVariables = controlledVariables;
    }

    public ExperimentalDesign(String name, String description, Set<Variable> responseVariables,
            Set<Variable> controlledVariables, Person author) {
        super(name, author, description);
        this.responseVariables = responseVariables;
        this.controlledVariables = controlledVariables;
    }

    public ExperimentalDesign(String name, String description, Set<Variable> responseVariables,
            Set<Variable> controlledVariables, BioDesign bioDesign, int site, String notes, Person author) {
        super(name, author, description);
        this.responseVariables = responseVariables;
        this.controlledVariables = controlledVariables;
        this.bioDesign = bioDesign;
        this.site = site;
        this.notes = notes;
    }

    public ExperimentalCondition createExperimentalCondition() {
        ExperimentalCondition experimentalCondition = new ExperimentalCondition();
        addExperimentalCondition(experimentalCondition);
        return experimentalCondition;
    }
    
    public void addExperimentalCondition(ExperimentalCondition experimentalCondition) {
    	if (experimentalConditions == null) {
    		experimentalConditions = new HashSet<ExperimentalCondition>();
    	}
    	experimentalConditions.add(experimentalCondition);
    }
    
    public void addSubDesign(ExperimentalDesign subDesign) {
    	if (subDesigns == null) {
    		subDesigns = new HashSet<ExperimentalDesign>();
    	}
    	subDesigns.add(subDesign);
    }
    
    public void addResponseVariable(Variable responseVariable) {
    	if (responseVariables == null) {
    		responseVariables = new HashSet<Variable>();
    	}
    	responseVariables.add(responseVariable);
    }
    
    public void addControlledVariable(Variable controlledVariable) {
    	if (controlledVariables == null) {
    		controlledVariables = new HashSet<Variable>();
    	}
    	controlledVariables.add(controlledVariable);
    }

    public Map getMap(){
        Map map = new HashMap();
        map.put("name", this.getName());
        map.put("author", this.getAuthor().getName());
        map.put("description", this.getDescription());
        if (this.responseVariables!=null)
            map.put("responseVariables", this.responseVariables.size());
        if (this.controlledVariables!=null)
            map.put("controlledVariables", this.controlledVariables.size());
        if (this.experimentalConditions!=null)
            map.put("experimentalConditions", this.experimentalConditions.size());
        if (this.bioDesign!=null)
            map.put("bioDesign", this.bioDesign.getName());
        if (this.subDesigns!=null)
            map.put("subDesigns", this.subDesigns.size());
        if (this.parentDesign!=null)
            map.put("parentDesign", this.parentDesign.getName());
        map.put("site", this.site);
        map.put("comments", this.notes);
        return map;
    }
    
    public JSONObject getJSON(){
        JSONObject obj = new JSONObject();
        obj.put("name", this.getName());
        obj.put("author", this.getAuthor().getName());
        obj.put("description", this.getDescription());
        if (this.responseVariables!=null)
            obj.put("responseVariables", this.responseVariables.size());
        if (this.controlledVariables!=null)
            obj.put("controlledVariables", this.controlledVariables.size());
        if (this.experimentalConditions!=null)
            obj.put("experimentalConditions", this.experimentalConditions.size());
        if (this.bioDesign!=null)
            obj.put("bioDesign", this.bioDesign.getName());
        if (this.subDesigns!=null)
            obj.put("subDesigns", this.subDesigns.size());
        if (this.parentDesign!=null)
            obj.put("parentDesign", this.parentDesign.getName());
        obj.put("site", this.site);
        obj.put("comments", this.notes);
        return obj;
    }
    
    public String toString(){
        String str = "";
        str += "Name : " + this.getName() + "\n";
        str += "Author : " + this.getAuthor().getName() + "\n";
        str += "Description : " + this.getDescription() + "\n";
        if (this.responseVariables!=null)
            str += "Response Variables : " + this.responseVariables.size() + "\n";
        if (this.controlledVariables!=null)
            str += "Controlled Variables : " + this.controlledVariables.size() + "\n";
        if (this.experimentalConditions!=null)
            str += "Experimental Conditions : " + this.experimentalConditions.size() + "\n";
        if (this.bioDesign!=null)
            str += "Bio Design : " + this.bioDesign.getName() + "\n";
        if (this.subDesigns!=null)
            str += "Sub Designs : " + this.subDesigns.size() + "\n";
        if (this.parentDesign!=null)
            str += "Parent Design : " + this.parentDesign.getName() + "\n";
        str += "Integration Site : " + this.site + "\n";
        str += "Additional Comments : " + this.notes;
        return str;
    }
}

package org.clothocad.model;

import java.util.HashMap;
import org.clothocad.core.persistence.annotations.ReferenceCollection;

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
public class BasicModule extends Module {

    @NotNull
    @Size(min=1)
    @Getter
    @Setter
    @ReferenceCollection
    protected Set<Feature> features;

    public BasicModule(String name, ModuleRole role, Set<Feature> features, Person author) {
        super(name, role, author);
        this.features = features;
    }

    public BasicModule(String name, String description, ModuleRole role, Set<Feature> features, Person author) {
        super(name, description, role, author);
        this.features = features;
    }
    
    public void addFeature(Feature feature) {
    	if (features == null) {
    		features = new HashSet<Feature>();
    	}
    	features.add(feature);
    }

    public Map getMap(){
        Map map = new HashMap();
        map.put("name", this.getName());
        map.put("author", this.getAuthor().getName());
        map.put("description", this.getDescription());
        if (features!=null) {
            int i = 0;
            for (Feature feature : features) {
                map.put("feature" + i, feature.getName());
                i++;
            }
        }
        map.put("role", this.role.toString());
        if (influences!=null) {
            int j = 0;
            for (Influence influence : influences) {
                map.put("influence" + j, influence.getName());
                j++;
            }
        }
        if (parentModule!=null) {
            map.put("name", this.parentModule.getName());
        }
        return map;
    }
    
    public JSONObject getJSON(){
        JSONObject obj = new JSONObject();
        obj.put("name", this.getName());
        obj.put("author", this.getAuthor().getName());
        obj.put("description", this.getDescription());
        if (features!=null) {
            int i = 0;
            for (Feature feature : features) {
                obj.put("feature" + i, feature.getName());
                i++;
            }
        }
        obj.put("role", this.role.toString());
        if (influences!=null) {
            int j = 0;
            for (Influence influence : influences) {
                obj.put("influence" + j, influence.getName());
                j++;
            }
        }
        if (parentModule!=null) {
            obj.put("name", this.parentModule.getName());
        }
        return obj;
    }
    
    public String toString(){
        String str = "";
        str += "Name : " + this.getName() + "\n";
        str += "Author : " + this.getAuthor().getName() + "\n";
        str += "Description : " + this.getDescription() + "\n";
        if (features!=null) {
            int i = 0;
            for (Feature feature : features) {
                str += "Feature " + i + " : " + feature.getName() + "\n";
                i++;
            }
        }
        str += "Module Role : " + this.role.toString();
        if (influences!=null) {
            int j = 0;
            for (Influence influence : influences) {
                str += "\nInfluence " + j + " : " + influence.getName();
                j++;
            }
        }
        if (parentModule!=null) {
            str += "\nParent Module : " + this.parentModule.getName();
        }
        return str;
    }
}

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.model;

import java.util.HashMap;
import java.util.Map;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.clothocad.core.datums.SharableObjBase;
import org.json.simple.JSONObject;

/**
 *
 * @author mardian
 */
@NoArgsConstructor
public class Expression extends SharableObjBase {
    
    @Setter
    @Getter
    protected double level;

    @Setter
    @Getter
    protected String software;

    @Setter
    @Getter
    protected Feature feature;

    @Getter
    @Setter
    protected BioDesign bioDesign;

    public Expression(String name, String description, Person author) {
        super(name, author, description);
    }
    
    public Expression(String name, String description, double level, String software, BioDesign bioDesign, Feature feature, Person author) {
        super(name, author, description);
        this.level = level;
        this.software = software;
        this.bioDesign = bioDesign;
        this.feature = feature;
    }
    
    public Map getMap(){
        Map map = new HashMap();
        map.put("name", this.getName());
        map.put("author", this.getAuthor().getName());
        map.put("description", this.getDescription());
        map.put("level", this.level);
        map.put("software", this.software);
        if(feature!=null)
            map.put("feature", this.feature.getName());
        if(bioDesign!=null)
            map.put("bioDesign", this.bioDesign.getName());
        return map;
    }
    
    public JSONObject getJSON(){
        JSONObject obj = new JSONObject();
        obj.put("name", this.getName());
        obj.put("author", this.getAuthor().getName());
        obj.put("description", this.getDescription());
        obj.put("level", this.level);
        obj.put("software", this.software);
        if(feature!=null)
            obj.put("feature", this.feature.getName());
        if(bioDesign!=null)
            obj.put("bioDesign", this.bioDesign.getName());
        return obj;
    }
    
    public String toString(){
        String str = "";
        str += "Name : " + this.getName() + "\n";
        str += "Author : " + this.getAuthor().getName() + "\n";
        str += "Description : " + this.getDescription() + "\n";
        str += "Level : " + this.level + "\n";
        str += "Software : " + this.software;
        if(feature!=null)
            str += "\nFeature" + this.feature.getName();
        if(bioDesign!=null)
            str += "\nBio Design" + this.bioDesign.getName();
        return str;
    }
}

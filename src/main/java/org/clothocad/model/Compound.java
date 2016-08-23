/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.model;

import java.util.HashMap;
import java.util.Map;
import javax.validation.constraints.NotNull;
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
public class Compound extends SharableObjBase {
    
    @NotNull
    @Setter
    @Getter
    protected CompoundType type;
    
    @Setter
    @Getter
    protected boolean isProduct;
    
    @Setter
    @Getter
    protected String pubChemId;
    
    @Setter
    @Getter
    protected double level;
    
    @Setter
    @Getter
    protected String formula;
    
    public static enum CompoundType {
    	METABOLITE, PROTEIN;
    }
    
    /**
     * Constructor of a new Compound
     *
     * @param name
     * @param description
     * @param author
     */
    public Compound(String name, String description, Person author) {
        super(name, author, description);
        this.type = CompoundType.PROTEIN;
    }
    
    /**
     * Constructor of a new Compound
     *
     * @param name
     * @param description
     * @param author
     * @param type
     */
    public Compound(String name, String description, Person author, CompoundType type) {
        super(name, author, description);
        this.type = type;
    }
    
    
    /**
     * Constructor of a new Compound
     *
     * @param name
     * @param description
     * @param author
     * @param type
     * @param isProduct
     * @param pubChemId
     * @param level
     * @param formula
     */
    public Compound(String name, String description, Person author, CompoundType type,
            boolean isProduct, String pubChemId, double level, String formula) {
        this(name, description, author, type);
        this.isProduct = isProduct;
        this.pubChemId = pubChemId;
        this.level = level;
        this.formula = formula;
    }
    
    public Map getMap(){
        Map map = new HashMap();
        map.put("name", this.getName());
        map.put("author", this.getAuthor().getName());
        map.put("description", this.getDescription());
        map.put("compoundType", this.type.toString());
        map.put("isProduct", this.isProduct);
        map.put("pubChemId", this.pubChemId);
        map.put("level", this.level);
        map.put("formula", this.formula);
        return map;
    }
    
    public JSONObject getJSON(){
        JSONObject obj = new JSONObject();
        obj.put("name", this.getName());
        obj.put("author", this.getAuthor().getName());
        obj.put("description", this.getDescription());
        obj.put("compoundType", this.type.toString());
        obj.put("isProduct", this.isProduct);
        obj.put("pubChemId", this.pubChemId);
        obj.put("level", this.level);
        obj.put("formula", this.formula);
        return obj;
    }
    
    public String toString(){
        String str = "";
        str += "Name : " + this.getName() + "\n";
        str += "Name : " + this.getName() + "\n";
        str += "Author : " + this.getAuthor().getName() + "\n";
        str += "Description : " + this.getDescription() + "\n";
        str += "Compound Type : " + this.type.toString() + "\n";
        if (this.isProduct)
            str += "Is a Product : Yes\n";
        else
            str += "Is a Product : No\n";
        str += "Identifier : " + this.pubChemId + "\n";
        str += "Level : " + this.level + "\n";
        str += "Formula : " + this.formula;
        return str;
    }
}

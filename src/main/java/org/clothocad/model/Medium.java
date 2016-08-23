package org.clothocad.model;

import java.util.HashMap;
import java.util.Map;
import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.Reference;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.json.simple.JSONObject;

/**
 * 
 * @author Nicholas Roehner
 * @author mardian
 */ 
@NoArgsConstructor
public class Medium extends SharableObjBase {

    /*@Getter
    @Setter
    @Reference
    protected String medium;*/

    @Getter
    @Setter
    @Reference
    protected Medium parentMedium;

    public Medium(String name, Person author) {
        super(name, author);
    }

    public Medium(String name, String description, Person author) {
        super(name, author, description);
    }

    /*public Medium(String name, String description, String medium, Person author) {
        super(name, author, description);
        this.medium = medium;
    }*/

    public Map getMap(){
        Map map = new HashMap();
        map.put("name", this.getName());
        map.put("author", this.getAuthor().getName());
        map.put("description", this.getDescription());
        //map.put("medium", this.medium);
        if(parentMedium!=null)
            map.put("parent", this.parentMedium.getName());
        return map;
    }
    
    public JSONObject getJSON(){
        JSONObject obj = new JSONObject();
        obj.put("name", this.getName());
        obj.put("author", this.getAuthor().getName());
        obj.put("description", this.getDescription());
        //obj.put("medium", this.medium);
        if(parentMedium!=null)
            obj.put("parent", this.parentMedium.getName());
        return obj;
    }
    
    public String toString(){
        String str = "";
        str += "Name : " + this.getName() + "\n";
        str += "Author : " + this.getAuthor().getName() + "\n";
        str += "Description : " + this.getDescription() + "\n";
        //str += "Medium : " + this.medium;
        if(parentMedium!=null)
            str +=  "\nParent medium : " + this.parentMedium.getName();
        return str;
    }
}

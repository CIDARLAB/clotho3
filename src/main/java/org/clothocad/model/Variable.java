package org.clothocad.model;

import java.util.HashMap;
import java.util.Map;
import lombok.NoArgsConstructor;

import org.clothocad.core.datums.SharableObjBase;
import org.json.simple.JSONObject;

/**
* @author Nicholas Roehner
* @author mardian
*/
@NoArgsConstructor
public class Variable extends SharableObjBase {

    public Variable(String name, Person author) {
        super(name, author);
    }

    public Variable(String name, String description, Person author) {
        super(name, author, description);
    }
    
    public Map getMap(){
        Map map = new HashMap();
        map.put("name", this.getName());
        map.put("author", this.getAuthor().getName());
        map.put("description", this.getDescription());
        return map;
    }
    
    public JSONObject getJSON(){
        JSONObject obj = new JSONObject();
        obj.put("name", this.getName());
        obj.put("author", this.getAuthor().getName());
        obj.put("description", this.getDescription());
        return obj;
    }
    
    public String toString(){
        String str = "";
        str += "Name : " + this.getName();
        str += "Author : " + this.getAuthor().getName();
        str += "Description : " + this.getDescription();
        return str;
    }

}

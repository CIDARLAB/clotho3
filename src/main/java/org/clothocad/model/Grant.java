/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.model;

import java.util.Date;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import lombok.Getter;
import lombok.Setter;
import org.clothocad.core.datums.SharableObjBase;
import org.json.simple.JSONObject;

/**
 *
 * @author mardian
 */
public class Grant extends SharableObjBase {
    
    @Setter
    @Getter
    private String grantName;
    
    @Setter
    @Getter
    private double amount;
    
    @Setter
    @Getter
    private Date start, finished;
    
    @Setter
    @Getter
    private Person PI;
    
    @Setter
    @Getter
    private Set<Lab> labs;
    
    @Setter
    @Getter
    private Institution sponsor;
    
    public Grant(String name, String description, Person author) {
        super(name, author, description);
    }
    
    public Grant(String name, String description, String sponsor, Person author) {
        super(name, author, description);
        this.sponsor = new Institution(sponsor, "", author);
    }
    
    public Map getMap(){
        Map map = new HashMap();
        //map.put("id", this.getId().getValue());
        map.put("name", this.getName());
        map.put("author", this.getAuthor().getName());
        map.put("description", this.getDescription());
        map.put("sponsor", this.sponsor.getName());
        return map;
    }
    
    public JSONObject getJSON(){
        JSONObject obj = new JSONObject();
        //obj.put("id", this.getId().getValue());
        obj.put("name", this.getName());
        obj.put("author", this.getAuthor().getName());
        obj.put("description", this.getDescription());
        obj.put("sponsor", this.sponsor.getName());
        return obj;
    }
    
    public String toString(){
        String str = "";
        //str += "ID : " + this.getId().getValue() + "\n";
        str += "Name : " + this.getName() + "\n";
        str += "Author : " + this.getAuthor().getName() + "\n";
        str += "Description : " + this.getDescription() + "\n";
        str += "Sponsor : " + this.sponsor.getName();
        return str;
    }
}

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
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.clothocad.core.datums.SharableObjBase;
import org.json.simple.JSONObject;

/**
 *
 * @author mardian
 */
@NoArgsConstructor
public class Publication extends SharableObjBase {
    
    @Setter
    @Getter
    private String doi, title, journal, publisher, url;
    
    @Setter
    @Getter
    private Project project;
    
    @Setter
    @Getter
    private Date submitted, accepted, published;
    
    @Setter
    @Getter
    private Person correspondingAuthor;
    
    @Setter
    @Getter
    private Set<Person> authors;
    
    @Setter
    @Getter
    private Set<Publication> references;
    
    @Setter
    @Getter
    private PublicationType type;
    
    @Setter
    @Getter
    private PublicationStatus status;
    
    public static enum PublicationType {
    	ARTICLE, REVIEW, NEWS, OPINION;
    }
    
    public static enum PublicationStatus {
    	IN_PREPARATION, SUBMITTED, UNDER_REVIEW, ACCEPTED, UNDER_REVISION, REJECTED, PUBLISHED;
    }
    
    public Publication(String name, String description, Person author) {
        super(name, author, description);
    }
    
    public Map getMap(){
        Map map = new HashMap();
        //map.put("id", this.getId().getValue());
        map.put("name", this.getName());
        map.put("author", this.getAuthor().getName());
        map.put("description", this.getDescription());
        return map;
    }
    
    public JSONObject getJSON(){
        JSONObject obj = new JSONObject();
        //obj.put("id", this.getId().getValue());
        obj.put("name", this.getName());
        obj.put("author", this.getAuthor().getName());
        obj.put("description", this.getDescription());
        return obj;
    }
    
    public String toString(){
        String str = "";
        //str += "ID : " + this.getId().getValue() + "\n";
        str += "Name : " + this.getName() + "\n";
        str += "Author : " + this.getAuthor().getName() + "\n";
        str += "Description : " + this.getDescription() + "\n";
        return str;
    }
}

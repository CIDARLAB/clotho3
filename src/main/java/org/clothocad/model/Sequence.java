package org.clothocad.model;

import java.util.HashMap;
import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.Reference;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import javax.validation.constraints.NotNull;
import javax.validation.constraints.Pattern;
import org.json.simple.JSONObject;

/**
* @author Nicholas Roehner
* @author mardian
*/
@NoArgsConstructor
public class Sequence extends SharableObjBase {

    @NotNull
    @Getter
    @Setter
    @Pattern(regexp="[ATUCGRYKMSWBDHVN]*", flags={Pattern.Flag.CASE_INSENSITIVE})
    protected String sequence;

    @Getter
    @Setter
    protected Set<Annotation> annotations;

    @Getter
    @Setter
    @Reference
    protected Sequence parentSequence;

    public Sequence(String name, String sequence, Person author) {
        super(name, author);
        this.sequence = sequence;
    }

    public Sequence(String name, String description, String sequence, Person author) {
        super(name, author, description);
        this.sequence = sequence;
    }

    public Annotation createAnnotation(String name, int start, int end, boolean isForwardStrand,
            Person author) {
        Annotation annotation = new Annotation(name, start, end, isForwardStrand, author);
        addAnnotation(annotation);
        return annotation;
    }

    public Annotation createAnnotation(String name, String description, int start, int end,
            boolean isForwardStrand, Person author) {
        Annotation annotation = new Annotation(name, description, start, end, isForwardStrand, author);
        addAnnotation(annotation);
        return annotation;
    }
    
    public void addAnnotation(Annotation annotation) {
    	if (annotations == null) {
    		annotations = new HashSet<Annotation>();
    	}
    	annotations.add(annotation);
    }

    public Map getMap(){
        Map map = new HashMap();
        map.put("name", this.getName());
        map.put("author", this.getAuthor().getName());
        map.put("description", this.getDescription());
        map.put("sequence", this.sequence);
        map.put("length", this.sequence.length());
        if(annotations!=null)
            map.put("annotation", this.annotations.size());
        if (parentSequence!=null)
            map.put("parent", this.parentSequence.getName());
        return map;
    }
    
    public JSONObject getJSON(){
        JSONObject obj = new JSONObject();
        obj.put("name", this.getName());
        obj.put("author", this.getAuthor().getName());
        obj.put("description", this.getDescription());
        obj.put("sequence", this.sequence);
        obj.put("length", this.sequence.length());
        if(annotations!=null)
            obj.put("annotation", this.annotations.size());
        if (parentSequence!=null)
            obj.put("parent", this.parentSequence.getName());
        return obj;
    }
    
    public String toString(){
        String str = "";
        str += "Name : " + this.getName() + "\n";
        str += "Author : " + this.getAuthor().getName() + "\n";
        str += "Description : " + this.getDescription() + "\n";
        str += "Sequence : " + this.sequence;
        if(annotations!=null)
            str +=  "\nNumber of annotations : " + this.annotations.size();
        if (parentSequence!=null)
            str += "\nParent sequence : " + this.parentSequence.getName();
        return str;
    }
}

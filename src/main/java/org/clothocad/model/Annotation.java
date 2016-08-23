package org.clothocad.model;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.Reference;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import java.awt.Color;
import java.util.HashMap;
import java.util.Map;

import javax.validation.constraints.NotNull;
import org.json.simple.JSONObject;

/**
 * An Annotation is a single line of Genbank essentially.  It maps a Feature
 * or just something a user has labeled as signficant on a Sequence.
 */

/**
 * 
 * @author J. Christopher Anderson
 * @author Nicholas Roehner
 * @author mardian
 */
@NoArgsConstructor
public class Annotation extends SharableObjBase {

    @Getter
    @Setter
    protected String symbol;

    @NotNull
    @Getter
    @Setter
    protected boolean isForwardStrand;

    @Getter
    @Setter
    @Reference
    protected Feature feature;

    @NotNull
    @Getter
    @Setter
    protected int start, end;

    @Setter
    protected Color forwardColor, reverseColor;

    /**
     * @param name
     * @param seq
     * @param end
     * @param start
     * @param author
     * @param isForwardStrand
     */
    public Annotation(String name, int start, int end, boolean isForwardStrand, Person author) {
        super(name, author);
        this.start = start;
        this.end = end;
        this.isForwardStrand = isForwardStrand;
    }

    /**
     * @param name
     * @param description
     * @param seq
     * @param end
     * @param start
     * @param author
     * @param isForwardStrand
     */
    public Annotation(String name, String description, int start, int end,
            boolean isForwardStrand, Person author) {
        super(name, author, description);
        this.start = start;
        this.end = end;
        this.isForwardStrand = isForwardStrand;
    }

    /**
     * Reverse the orientation of the annotation (reverse complement
     * it and flip flop the start and end sites).  Called from NucSeq
     * when it's reverse complemented
     * @param seqLength
     */
    public void invert(int seqLength) {
        isForwardStrand = !isForwardStrand;
        int s = start;
        start = seqLength - end;
        end = seqLength - s;
    }

    /**
     * Get the approriate color for the annoation
     * @return either the forward or reverse color depending
     * on the orientation of the annotation
     */
    public Color getColor() {
        if (isForwardStrand) {
            return getForwardColor();
        } else {
            return getReverseColor();
        }
    }

    /**
     * Get the forward color as an integer code
     * @return an integer of the Color
     */
    public int getForwardColorAsInt() {
        return getForwardColor().getRGB();
    }

    /**
     * Get the reverse color as an integer code
     * @return an integer of the Color
     */
    public int getReverseColorAsInt() {
        return getReverseColor().getRGB();
    }

    /**
     * Get the preferred forward color for this Annotation.  If no forward color
     * was set, a default color will be returned.
     * @return an AWT Color object.  It won't be null;
     */
    public Color getForwardColor() {
        if (forwardColor == null) {
            forwardColor = new Color(125, 225, 235);
        }
        return forwardColor;
    }

    /**
     * Get the preferred reverse color for this Annotation.  If no reverse color
     * was set, a default color will be returned.
     * @return an AWT Color object.  It won't be null;
     */
    public Color getReverseColor() {
        if (reverseColor == null) {
            reverseColor = new Color(125, 225, 235);
        }
        return reverseColor;
    }

    /**
     * Set the forward and reverse preferred colors for this feature to some
     * random medium-bright color.
     */
    public void setRandomColors() {
        int[][] intVal = new int[2][3];
        for (int j = 0; j < 3; j++) {
            double doubVal = Math.floor(Math.random() * 155 + 100);
            intVal[0][j] = (int) doubVal;
            intVal[1][j] = 255 - intVal[0][j];
        }
        forwardColor = new Color(intVal[0][0], intVal[0][1], intVal[0][2]);
        reverseColor = new Color(intVal[1][0], intVal[1][1], intVal[1][2]);
    }

    public Map getMap(){
        Map map = new HashMap();
        map.put("name", this.getName());
        map.put("author", this.getAuthor().getName());
        map.put("description", this.getDescription());
        map.put("symbol", this.symbol);
        map.put("isForwardStrand", this.isForwardStrand);
        if (feature!=null)
            map.put("feature", this.feature.getName());
        map.put("start", this.start);
        map.put("end", this.end);
        if (forwardColor!=null)
            map.put("forwardColor", this.forwardColor.toString());
        if (reverseColor!=null)
            map.put("reverseColor", this.reverseColor.toString());
        return map;
    }
    
    public JSONObject getJSON(){
        JSONObject obj = new JSONObject();
        obj.put("name", this.getName());
        obj.put("author", this.getAuthor().getName());
        obj.put("description", this.getDescription());
        obj.put("symbol", this.symbol);
        obj.put("isForwardStrand", this.isForwardStrand);
        if (feature!=null)
            obj.put("feature", this.feature.getName());
        obj.put("start", this.start);
        obj.put("end", this.end);
        if (forwardColor!=null)
            obj.put("forwardColor", this.forwardColor.toString());
        if (reverseColor!=null)
            obj.put("reverseColor", this.reverseColor.toString());
        return obj;
    }
    
    public String toString(){
        String str = "";
        str += "name" + this.getName() + "\n";
        str += "author" + this.getAuthor().getName() + "\n";
        str += "description" + this.getDescription() + "\n";
        str += "symbol" + this.symbol + "\n";
        str += "isForwardStrand" + this.isForwardStrand + "\n";
        if (feature!=null)
            str += "feature" + this.feature.getName() + "\n";
        str += "start" + this.start + "\n";
        str += "end" + this.end;
        if (forwardColor!=null)
            str += "\nforwardColor" + this.forwardColor.toString();
        if (reverseColor!=null)
            str += "\nreverseColor" + this.reverseColor.toString();
        return str;
    }
}

package org.clothocad.model;

import java.util.HashMap;
import java.util.Map;
import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.Reference;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import javax.validation.constraints.NotNull;
import org.json.simple.JSONObject;

/**
 *
 * @author J. Christopher Anderson
 * @author Nicholas Roehner
 * @author mardian
 */
@NoArgsConstructor
public class Feature extends SharableObjBase {

    @Setter
    @Getter
    @Reference
    protected Sequence sequence;

    @Setter
    @Getter
    protected String genbankId, swissProtId, PDBId;

    @Getter
    protected short riskGroup;

    @NotNull
    @Setter
    @Getter
    protected FeatureRole role;

    @Getter
    @Setter
    @Reference
    protected Feature parentFeature;

    /**
     * Constructor of a new Feature
     *
     * @param name
     * @param role
     * @param author
     */
    public Feature(String name, FeatureRole role, Person author) {
        super(name, author);
        this.role = role;
    }

    /**
     * Constructor of a new Feature
     *
     * @param name
     * @param description
     * @param sequence
     * @param role
     * @param author
     */
    public Feature(String name, String description, Sequence sequence, FeatureRole role, Person author) {
        super(name, author, description);
        this.sequence = sequence;
        this.role = role;
    }

    /**
     * Constructor of a new Feature
     *
     * @param name
     * @param description
     * @param role
     * @param author
     */
    public Feature(String name, String description, FeatureRole role, Person author) {
        super(name, author, description);
        this.role = role;
    }

    /**
     * Change the risk group of the Feature. You can only raise the risk group.
     *
     * @param newrg the new risk group (1 through 5)
     */
    public final void setRiskGroup(short newrg) {
        if (newrg > riskGroup && newrg <= 5) {
            //addUndo("_riskGroup", _featDatum._riskGroup, newrg);
            riskGroup = newrg;
            // setChanged(org.clothocore.api.dnd.RefreshEvent.Condition.RISK_GROUP_CHANGED);
        }
        //todo: throw appropriate invalid operation exception
    }

    // Feel free to add more of these
    public static enum FeatureRole {
    	BARCODE, CDS, DEGRADATION_TAG, GENE, LOCALIZATION_TAG, OPERATOR, PROMOTER, SCAR, SPACER, RBS, RIBOZYME, TERMINATOR,
        TOXICITY_TEST;
    }
    
    public Map getMap(){
        Map map = new HashMap();
        map.put("name", this.getName());
        map.put("author", this.getAuthor().getName());
        map.put("description", this.getDescription());
        map.put("sequence", this.sequence.getSequence());
        map.put("genbankId", this.genbankId);
        map.put("swissProtId", this.swissProtId);
        map.put("PDBId", this.PDBId);
        map.put("riskGroup", this.riskGroup);
        map.put("role", this.role.toString());
        if (this.parentFeature!=null)
            map.put("parent", this.parentFeature.getName());
        return map;
    }
    
    public JSONObject getJSON(){
        JSONObject obj = new JSONObject();
        obj.put("name", this.getName());
        obj.put("author", this.getAuthor().getName());
        obj.put("description", this.getDescription());
        obj.put("sequence", this.sequence.getSequence());
        obj.put("genbankId", this.genbankId);
        obj.put("swissProtId", this.swissProtId);
        obj.put("PDBId", this.PDBId);
        obj.put("riskGroup", this.riskGroup);
        obj.put("role", this.role.toString());
        if (this.parentFeature!=null)
            obj.put("parent", this.parentFeature.getName());
        return obj;
    }
    
    public String toString(){
        String str = "";
        str += "Name : " + this.getName() + "\n";
        str += "Author : " + this.getAuthor().getName() + "\n";
        str += "Description : " + this.getDescription() + "\n";
        str += "Sequence : " + this.sequence.getSequence() + "\n";
        str += "GenBank ID : " + this.genbankId + "\n";
        str += "SwissProt ID : " + this.swissProtId + "\n";
        str += "PDB ID : " + this.PDBId + "\n";
        str += "Risk Group : " + this.riskGroup + "\n";
        str += "Feature Role : " + this.role.toString();
        if (this.parentFeature!=null)
            str += "\nParent Feature : " + this.parentFeature.getName();
        return str;
    }
}

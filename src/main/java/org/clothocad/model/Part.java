/*
Copyright (c) 2009 The Regents of the University of California.
All rights reserved.
Permission is hereby granted, without written agreement and without
license or royalty fees, to use, copy, modify, and distribute this
software and its documentation for any purpose, provided that the above
copyright notice and the following two paragraphs appear in all copies
of this software.

IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.
createQuery
THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
ENHANCEMENTS, OR MODIFICATIONS..
 */
package org.clothocad.model;

import com.github.jmkgreen.morphia.annotations.Reference;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import javax.validation.Valid;
import javax.validation.constraints.AssertTrue;
import javax.validation.constraints.NotNull;

import org.clothocad.core.datums.ObjBase;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import lombok.ToString;
import lombok.extern.slf4j.Slf4j;
import org.clothocad.core.persistence.Replace;

/**
 *
 * @author jcanderson
 */
@ToString(callSuper=true, includeFieldNames=true)
@NoArgsConstructor
@Slf4j
public class Part extends ObjBase {
    
    //inject connection

    
    @Setter
    @Getter
    @Reference
    private Person author;
    
    @Getter
    @NotNull
    @Replace(encoder = "getFormatName", decoder="setFormatFromName")
    private Format format;

    @Getter
    @Setter
    private String shortDescription;
    
    @Valid
    @Reference
    private NucSeq sequence;
    
    @Getter
    @NotNull
    private PartType partType;
    
    @Getter
    @Setter
    @Reference
    private List<Part> composition;
    
    @Getter
    private short riskGroup;
    
    private Part(String name, String desc, Format form, Person author, PartType type){
        super(name); 
        this.shortDescription = desc;
        this.partType = type;
        this.author = author;
        this.format = form;

    }
        
     /**
     * Create basic from scratch
     *
     * @param name  nickname of Part, like "roo40"
     * @param shortdescription short description, such as "[TetR]"
     * @param seq sequence of the Part like "cgaaggcaggacacacatg"
     * @param form Part Format
     * @param author author of the Part
     * @param partType Basic or Composite
     */
    protected Part(String name, String shortDescription, String seq, Format form, Person author) {
        this(name, shortDescription, form, author, PartType.BASIC);
        this.sequence = new NucSeq(seq);
    }

    /**
     * Create composite from scratch
     *
     * @param name  nickname of Part, like "roo40"
     * @param shortdescription short description, such as "[TetR]"
     * @param seq sequence of the Part like "cgaaggcaggacacacatg"
     * @param form Part Format
     * @param author author of the Part
     * @param partType Basic or Composite
     */
    protected Part(String name, String shortdescription, Format form, Person author) {
        //Use the transient constructor and add to Collector later if the Part passes checks
        this(name, shortdescription, form, author, PartType.COMPOSITE);
    }

    /**
     * Call this method to construct a new basic Part.  It will check that:
     *  1) A sequence in this Format isn't already in the database
     *      otherwise it returns the Part that already exists
     *  2) The Part obeys its Format standard
     *      otherwise it returns null
     *
     * Only if those are satisfied is it added to the Collector and returned to
     * the calling code.
     *
     * @param name  nickname of Part, like "roo40"
     * @param shortdescription short description, such as "[TetR]"
     * @param seq sequence of the Part like "cgaaggcaggacacacatg"
     * @param form Part Format
     * @param author author of the Part
     * @param partType Basic or Composite
     */
    
    
    //TODO: fix validator classpath issues
    /*{{
      validator = Validation.buildDefaultValidatorFactory().getValidator();
    }}

    private static Validator validator;
    */
    
    //create part, validate, return null if bad
    //default shortdesc for composites is concatenation of component parts' shortdesc
    //copy composition array
    public static Part generateBasic(String name, String shortdescription, String seq, Format form, Person author) {
        Part part = new Part(name, shortdescription, seq, form, author);
        /*Set<ConstraintViolation<Part>> violations = validator.validate(part);
        
        if (violations.isEmpty()) return part;
        return null;*/
        return part;
    };
    
    
    //additional requirements are defined by particular format
    public static Part generateComposite(List<Part> composition, Object additionalRequirements, Format f, Person author, String name, String shortdescription) {
        if (!f.checkComposite(composition, additionalRequirements)) {
            System.out.println("generateComposite: Doesn't obey format, return null");
            return null;
        }
        
        Part part = new Part(name, shortdescription, f, author);
        part.setComposition(composition);
        /*Set<ConstraintViolation<Part>> violations = validator.validate(part);
        
        if (violations.isEmpty()) return part;
        return null;*/
        return part;
     
    }

    
    @AssertTrue
    public boolean checkFormat() {
        if (partType == PartType.BASIC){
            return format.checkPart(this);
        }
        else {
            return format.checkComposite(this.composition, null);
        }
    }
    
    public String getFormatName(){
        return format.getClass().getSimpleName();
    }

    public void setFormatFromName(Map<String, Object> dbObject){
        String name = (String) dbObject.get("format");
        try {
            //XXX: stupid hack
            //TODO: search schemas instead
            format = (Format) Class.forName("org.clothocad.model." + name).newInstance();
            
        } catch (ClassNotFoundException | InstantiationException | IllegalAccessException ex) {
            log.error("Couldn't find format {}", "org.clothocad.model." + name);
        } 
    }
    
    /* SETTERS
     * */

    /**
     * This is a convenience method, the real change to the sequence
     * happens in the linked NucSeq
     * @param newseq
     */
    public void setSequence(final String newseq) {
        if (newseq.equals("") || newseq == null) {
            return;
        }

        if (partType.equals(PartType.COMPOSITE)) {
            return;
        }

        final String oldseq = sequence.toString();

        sequence.APIchangeSeq(newseq);

        boolean isok = format.checkPart(this);
        if (!isok) {
            sequence.APIchangeSeq(oldseq);
            return;
        }

        //Change the risk group
        //riskGroup = sequence.performBiosafetyCheck();
    }

    /**
     * Change the Format of the Part
     * @param f
     */
    public void setFormat(Format f) {
        if (f == null) {
            return;
        }

        boolean ok = f.checkPart(this);
        if (!ok) {
            return;
        }

        format = f;
    }


    /* GETTERS
     * */
    /*public List<ObjBase> getPlasmids() {
        ClothoConnection c = Collector.getDefaultConnection();
        ClothoQuery mainQuery = c.createQuery(ObjType.PLASMID);
        ClothoQuery partQuery = mainQuery.createAssociationQuery(Plasmid.Fields.PART);
        partQuery.add(partQuery.getMatchesCrit(Part.Fields.NAME, this.getName()));
        List<ObjBase> results = mainQuery.getResults();
        return results;
    }*/


    public NucSeq getSequence() {
        if (partType.equals(PartType.BASIC)) {
            return sequence;
        } else {
            //cache seq?
            return format.generateCompositeSequence(composition, null);
        }
    }

    /*public final void changeRiskGroup(Short newrg) {
        if (newrg > _partDatum._riskGroup) {
            addUndo("_riskGroup", _partDatum._riskGroup, newrg);
            _partDatum._riskGroup = newrg;
        }
        setChanged(RefreshEvent.Condition.RISK_GROUP_CHANGED);
    }

    /*
     * Determines the risk group of the Part,
     * relayed from the initial call to NucSeq's
     * call to foreign server
     * @param rg
     *
    private void setRiskGroup(short rg) {
        if (rg == 5) {
            relayRiskGroup((short) 5);
            return;
        }
        if (this._partDatum._partType.equals(partType.Composite)) {
            short currentHighest = rg;
            boolean firsthigher = false;
            for (Part p : this.getCompositeParts()) {
                //If the Part's risk group hasn't been determined, this one isn't either
                if (p.getRiskGroup() == -1) {
                    relayRiskGroup((short) -1);
                    return;
                }

                //If a subpart has a 2+ risk group, increment highest
                if (p.getRiskGroup() > currentHighest) {
                    currentHighest = p.getRiskGroup();
                    relayRiskGroup(currentHighest);
                }

                
                //If a subpart has a 2+ risk group
                if (p.getRiskGroup() > 1) {
                    if (firsthigher) {
                        
                        //Throw a dialog asking for user to put in the new risk group
                        ButtonGroup group = new javax.swing.ButtonGroup();
                        String msgString = "This composite part joins two subparts with risk groups of 2 or higher.  What should the new value be?";
                        int numelements = 5 - currentHighest;
                        Object[] array = new Object[numelements + 1];
                        JRadioButton[] buttons = new JRadioButton[numelements];
                        for (short i = 0; i < numelements; i++) {
                            buttons[i] = new javax.swing.JRadioButton("Risk Group " + (i + currentHighest));
                            group.add(buttons[i]);
                            array[i + 1] = buttons[i];
                        }
                        array[0] = msgString;

                        /* // sbhatia commented this out
                        int sel = -1;
                        while (sel != 0) {
                            sel = JOptionPane.showConfirmDialog(null, array, "", JOptionPane.DEFAULT_OPTION, JOptionPane.QUESTION_MESSAGE);
                        }*/
                        /*buttons[0].setSelected(true); // sbhatia added this line
                        scanButtons:
                        for (short i = 1; i < numelements; i++) {
                            if (buttons[i].isSelected()) {
                                relayRiskGroup((short) (currentHighest + 1));
                                break scanButtons;
                            }
                        }
                        currentHighest = _partDatum._riskGroup;

                    } else {
                        firsthigher = true;
                    }
                }
            }

            //If it's a basic Part, the risk group is whatever the algorithm said
        } else {
            relayRiskGroup(rg);
        }

        //Check it's features to see if any are higher RG
        for (Annotation an : this.getSeq().getAnnotations()) {
            Feature afeat = an.getFeature();
            if (afeat == null) {
                continue;
            }
            relayRiskGroup(afeat.getRiskGroup());
        }
    }

    private void relayRiskGroup(short value) {
        if (value > _partDatum._riskGroup) {
            _partDatum._riskGroup = value;
            System.out.println("Setting risk group to " + _partDatum._riskGroup);
            setChanged(RefreshEvent.Condition.RISK_GROUP_CHANGED);
        } else {
            fireData(new RefreshEvent(this, RefreshEvent.Condition.RISK_GROUP_CHANGED));
        }
    }*/

    public static Part retrieveByName(String name) {
        //query connection for one part whose name contains the provided string    
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public static Part retrieveByExactName(String name) {
        //query connection for one part whose name is the provided string
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @AssertTrue
    private boolean checkDBConstraints(){
        //name must be unique
        //format(by name) + sequence(by content)combination should be unique
        return true;
    }

    public static enum PartType {
        BASIC, COMPOSITE
    };
}

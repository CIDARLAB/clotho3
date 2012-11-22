/*
 * 
Copyright (c) 2010 The Regents of the University of California.
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

THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
ENHANCEMENTS, OR MODIFICATIONS..
 */

package org.clothocad.core.aspects.Proctor;

import flexjson.JSONSerializer;
import java.util.ArrayList;
import java.util.List;
import org.clothocad.core.aspects.Collector;
import org.clothocad.core.datums.Sharable;
import org.clothocad.core.datums.objbases.Person;
import org.clothocad.core.datums.util.ClothoDate;
import org.clothocad.core.datums.util.Permissions;
import org.json.JSONObject;


/**
 * @author John Christopher Anderson
 */
public abstract class Paver implements Sharable {

    public abstract JSONObject makeCommandList() throws Exception;
    public abstract JSONObject makeTocJSON() throws Exception;
    
    @Override
    public JSONObject toJSON() {
        try {
            JSONSerializer serializer = new JSONSerializer().exclude("*.class");
            serializer.prettyPrint(true);
            String serial = serializer.deepSerialize( this );
            return new JSONObject(serial);
        } catch (Exception ex) {
            return null;
        }
    }

    @Override
    public Person extractAuthor() {
        Person out = (Person) Collector.get().getDatum(authorId);
        return out;
    }

    @Override
    public String getId() {
        return id;
    }
    
    /**
     * For the Yay's, Nay's, and usage of this Paver, calculate
     * an overall rating on the 0..5 scale
     * @return 
     */
    protected final double calculateRating() {
        int totalVotes = numYays + numNays;
        if(totalVotes==0) {
            return -1;
        }
        return numYays/totalVotes;
    }
    
    /**
     * Returns the typical duration in seconds
     * @return 
     */
    protected final int calculateDuration() {
        long[] durations = new long[sessionRecords.size()];
        for(int i=0; i<sessionRecords.size(); i++) {
            SessionRecord record = sessionRecords.get(i);
            durations[i] = record.timeFinished.getAbsolute() - record.timeInitiated.getAbsolute();
        }
        
        //Ideally this would do something fancier, but I have it doing the average
        long sum = 0;
        for(long value : durations) {
            sum+=value;
        }
        double daverage = Math.floor(0.0001 * (sum/durations.length)); //value in seconds
        return (int) daverage;
    }

    //Permissions
    private Permissions permissions = new Permissions();
    
    //Metadata
    protected String id;
    protected String title;
    protected String description;
    protected String authorId;
    protected String smallIconURL;
    protected String largeIconURL;
    
    private int numYays = 0;
    private int numNays = 0;
    private int timesAccessed = 0;
    private List<SessionRecord> sessionRecords = new ArrayList<SessionRecord>();
    
    private ClothoDate dateCreated = new ClothoDate();
    private ClothoDate dateLastModified = new ClothoDate();
    private ClothoDate dateLastAccessed = new ClothoDate();
}

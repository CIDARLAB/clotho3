/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.aspects.Proctor;

<<<<<<< HEAD
import java.util.ArrayList;
import java.util.List;
import org.json.JSONObject;

=======
import lombok.AllArgsConstructor;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
>>>>>>> f3b08d8fff3a3286cafa81e37305bac44897688f

/**
 *
 * @author jcanderson
 */
<<<<<<< HEAD
public abstract class Paver  {
    
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
            durations[i] = record.timeFinished.getTime()- record.timeInitiated.getTime();
        }
        
        //Ideally this would do something fancier, but I have it doing the average
        long sum = 0;
        for(long value : durations) {
            sum+=value;
        }
        double daverage = Math.floor(0.0001 * (sum/durations.length)); //value in seconds
        return (int) daverage;
    }

    public abstract JSONObject getContent();
    
    //Metadata
    protected String paver_title;
    protected PaverType type;
    
    
    private int numYays = 0;
    private int numNays = 0;
    private int timesAccessed = 0;
    private List<SessionRecord> sessionRecords = new ArrayList<SessionRecord>();
    
    public static enum PaverType {template, quiz}
=======
@AllArgsConstructor
@NoArgsConstructor
public abstract class Paver {
    @Getter
    @Setter
    private String paver_title;
    @Getter
    private String type;
>>>>>>> f3b08d8fff3a3286cafa81e37305bac44897688f
}


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
import org.clothocad.core.aspects.Collector;
import org.clothocad.core.datums.Doo;
import org.clothocad.core.datums.Sharable.SharableType;
import org.clothocad.core.datums.objbases.Person;
import org.clothocad.core.datums.util.ClothoDate;
import org.clothocad.core.datums.util.Permissions;
import org.clothocad.core.util.Logger;
import org.json.JSONArray;
import org.json.JSONObject;


/**
 * @author John Christopher Anderson
 */


public class Trail 
		extends Paver {

    public Trail(Person author) {
		super(author, SharableType.TRAIL);
	}

	/**
     * This is the commandList for the entire trail, which basically means relay
     * the first Paver's makeCommandList to display the first paver
     * @return 
     * @author John Christopher Anderson
     * @status correct as of 10/11/12
     */
    @Override
    public JSONObject makeCommandList() throws Exception {
        if(pavers.isEmpty()) {
            return null;
        }
        
        Paver first = fetchPaver(0);
        return first.makeCommandList();        
    }
    
    /**
     * Helper method for fetching Paver objects
     * @param index
     * @return 
     * @author John Christopher Anderson
     * @status correct as of 10/11/12
     */
    private Paver fetchPaver(int index) {
        String uuid = pavers.get(index);
        Paver first = (Paver) Collector.get().getDatum(uuid);
        return first;
    }

    /**
     * This should be a recursive call to all the Pavers gathering
     * up their JSON and returning it as a bundled up object
     * 
     * @return 
     * @author John Christopher Anderson
     * @status partially written
     */
    @Override
    public JSONObject makeTocJSON() throws Exception {
       JSONObject displayData = new JSONObject();
       displayData.put("title", ""); //TODO:  FIGURE OUT WHERE TO GET THIS INFO CORRECTLY
       displayData.put("description", ""); //TODO:  FIGURE OUT WHERE TO GET THIS INFO CORRECTLY
       //PROBABLY A FEW MORE METADATA FIELDS
       
       
       //Bundle up all the Paver's TOC data
       JSONArray tocArray = new JSONArray();
       for(int i=0; i<pavers.size(); i++) {
           Paver paver = fetchPaver(i);
           if(paver==null) {
               Logger.log(Logger.Level.WARN, "Error fetching a");
               continue;
           }
           tocArray.put(paver.makeTocJSON());
       }
       displayData.put("toc_array", tocArray);
       return displayData;
    }

    //Trail-specific fields
    private ArrayList<String> pavers = new ArrayList<String>();  //Is really a List<Paver>, loose-coupled
}

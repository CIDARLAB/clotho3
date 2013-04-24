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

import java.util.ArrayList;
import org.clothocad.model.Person;

import org.json.JSONObject;


/**
 * @author John Christopher Anderson
 */


public class Content 
		extends Paver {
	

    public Content(Person author) {
		super(author);
	}

	@Override
    public JSONObject makeCommandList() throws Exception {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public JSONObject makeTocJSON() throws Exception {
       JSONObject displayData = new JSONObject();
       displayData.put("title", this.title);
       displayData.put("description", this.description);
       displayData.put("author_id", this.getAuthor().getUUID());
       displayData.put("rating", this.calculateRating());
       //PROBABLY A FEW MORE METADATA FIELDS
       
       return displayData;
    }
    
    public class ViewElements {
        public String viewId;
        public String instanceId;
    }
    
    //Content-specific fields
    private ArrayList<ViewElements> viewElements = 
    		new ArrayList<ViewElements>();  //Pairings of Views and Instance data to stuff in them
}

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
import java.util.List;
import java.util.UUID;
import org.clothocad.core.datums.Doo;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.View;
import org.clothocad.model.Person;
import org.json.JSONArray;
import org.json.JSONObject;

/**
 * A quiz is a Sharable object in the sense that passing the data around
 * is allowed, and you can do CRUD operations on Quizes from the web
 * 
 * However, 'execution' of the quiz, or taking the Quiz, happens via the
 * Proctor and managed carefully.  It's not just passed to the client and
 * done in js.  It's done on the server, but the data itself is sometimes
 * going to be passed over the web for people to author these.
 * 
 * @author John Christopher Anderson
 */


public class Quiz 
		extends Paver {

    private Quiz(Person author,
    		String name, 
    		String description,
    		View view,
    		List<Question> questions,
    		boolean doRandom,
    		Function rubric) {
		super(author);

        this.description = description;
        this.viewId = view.getUUID().toString();
        this.questions = questions;
        this.doRandom = doRandom;
        this.rubric = rubric;
	}


	public static Quiz create(Person author,
                                String name, 
                                String description,
                                View view,
                                List<Question> questions,
                                boolean doRandom,
                                Function rubric) {

        try {
            //CHECK THE DATA FOR WELL-FORMEDNESS

            //Create the Quiz object
            Quiz quiz = new Quiz(author, name, description, view, questions, doRandom, rubric);
            return quiz;
        } catch(Exception err) {
            err.printStackTrace();
            return null;
        }
    }


    public boolean isDoRandom() {
        return doRandom;
    }

    public List<Question> getQuestions() {
        return questions;
    }

    public String getSmallIconURL() {
        return smallIconURL;
    }

    public String getViewId() {
        return viewId;
    }

    public String getBadgeId() {
        return badgeId;
    }

    public Function getRubric() {
        return rubric;
    }  
    
    @Override
    public JSONObject makeCommandList() throws Exception {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public JSONObject makeTocJSON() throws Exception {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    public static class Question {
        public JSONObject question;  //A List of these objects needs to be understandable to the view with viewId
        public JSONObject answer;    //Structure of the correct answer, or whatever the rubric needs to know to grade a submission
    }
    
    //Quiz-specific fields
    private String viewId;     //The view to show the content of the quiz
    private List<Question> questions;  //the data and answers that go into the quiz
    private boolean doRandom;  //Should select the questions randomly?
    private String badgeId;    //The badge that can be earned by passing this quiz
    private Function rubric;   //This reference to a Function needs to input the data from the client and output some score
    
}

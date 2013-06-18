/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.aspects.Proctor;

import org.clothocad.core.datums.Doo;
import org.clothocad.core.layers.communication.mind.Mind;
import org.json.JSONObject;

/**
 *
 * @author jcanderson
 */
public class QuizDoo extends Doo  {
    
    /**
     * The quizDoo is constructed from a quizJSON that describes the questions
     * and answers.  The Quiz object holds the state of what question/answer combo
     * has been send to the user.  The Quiz is more like a Session object holding state.
     * 
     *
     * @param quizJSON
     * @param mind 
     */
    public QuizDoo(JSONObject quizJSON, Mind mind, Doo parent) {
        super(parent, false);
        
    }
}

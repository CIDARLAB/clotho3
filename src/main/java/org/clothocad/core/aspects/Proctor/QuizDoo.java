/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.aspects.Proctor;

import java.util.Map;
import org.clothocad.core.datums.Doo;
import org.clothocad.core.execution.Mind;

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
    public QuizDoo(Map<String, Object> quizJSON, Mind mind, Doo parent) {
        super(parent, false);
        
    }
}

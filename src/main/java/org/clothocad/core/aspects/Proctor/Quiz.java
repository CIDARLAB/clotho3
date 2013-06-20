/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.aspects.Proctor;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

/**
 *
 * @author jcanderson
 */
@NoArgsConstructor
public abstract class Quiz extends Paver {

    public Quiz(String paver_title, String title,  String quiz_type, String question) {
        super(paver_title, "quiz");
        this.title = title;
        this.type = quiz_type;
        this.question = question;
    }

    @Getter
    @Setter
    private String title;;
    @Getter
    @Setter
    private String type;
    @Getter
    @Setter
    private String question;
}

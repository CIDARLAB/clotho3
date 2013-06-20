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
public class TemplatePaver extends Paver {

    public TemplatePaver(String paver_title, String template) {
        super(paver_title, "template");
        this.template = template;
    }

    @Getter
    @Setter
    private String template;
}

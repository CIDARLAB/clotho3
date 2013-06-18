/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.aspects.Proctor;

import lombok.Getter;
import lombok.Setter;

/**
 *
 * @author jcanderson
 */
public class TemplatePaver extends Paver {

    public TemplatePaver(String paver_title, String template) {
        super(paver_title, "template");
    }

    @Getter
    @Setter
    private String template;
}

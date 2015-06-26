/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.persistence;

import lombok.Getter;
import lombok.Setter;

import org.clothocad.model.Part;

/**
 *
 * @author spaige
 */
public class ExtendedPart extends Part {
    @Getter
    @Setter
    private String additionalParameters;
    
}

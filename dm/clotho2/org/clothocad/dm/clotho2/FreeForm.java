/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.dm.clotho2;

import java.util.List;

/**
 *
 * @author spaige
 */
public class FreeForm implements Format {

    @Override
    public boolean checkPart(Part p) {
        return true;
    }

    @Override
    public boolean checkComposite(List<Part> composition, Object additionalRequirements) {
        return true;
    }

    @Override
    public NucSeq generateCompositeSequence(List<Part> composition, Object additionalRequirements) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
}

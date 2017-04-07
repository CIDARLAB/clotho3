package org.clothocad.model;

import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.persistence.annotations.Reference;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import javax.validation.constraints.NotNull;

/**
*
* @author Nicholas Roehner
*/
@NoArgsConstructor
public class Level extends ObjBase {

    @NotNull
    @Getter
    protected Parameter parameter;

    @Getter
    @Setter
    @Reference
    protected BioDesign bioDesign;

    @Getter
    @Setter
    protected String description;

    protected Level(String name) {
        super(name);
    }

    public Parameter createParameter(String name, double value, String variable) {
        parameter = new Parameter(name, value, variable);
        return parameter;
    }

}

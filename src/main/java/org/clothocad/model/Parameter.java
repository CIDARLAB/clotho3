package org.clothocad.model;

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
public class Parameter {

    @NotNull
    @Getter
    @Setter
    protected double value;
    
    @NotNull
    @Getter
    @Setter
    protected String name;

    @NotNull
    @Getter
    @Setter
    @Reference
    protected String variable;

//    @NotNull
    @Getter
    @Setter
    @Reference
    protected String units;

    @Getter
    @Setter
    @Reference
    protected Derivation derivation;

    protected Parameter(String name, double value, String variable) {
        this.name = name;
        this.value = value;
        this.variable = variable;
    }

}

package org.clothocad.model;

import lombok.Getter;
import lombok.NoArgsConstructor;

import java.util.HashSet;
import java.util.Set;

import javax.validation.constraints.NotNull;
import javax.validation.constraints.Size;

/**
*
* @author Nicholas Roehner
*/
@NoArgsConstructor
public class ExperimentalCondition {

    @NotNull
    @Size(min=1)
    @Getter
    protected Set<Parameter> parameters;

    public Parameter createParameter(String name, double value, String variable) {
        Parameter parameter = new Parameter(name, value, variable);
        addParameter(parameter);
        return parameter;
    }

    public void addParameter(Parameter parameter) {
    	if (parameters == null) {
    		parameters = new HashSet<Parameter>();
    	}
    	parameters.add(parameter);
    }
}

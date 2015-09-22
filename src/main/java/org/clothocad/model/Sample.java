package org.clothocad.model;

import java.util.HashSet;
import java.util.Set;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.Reference;
import org.clothocad.core.persistence.annotations.ReferenceCollection;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import javax.validation.constraints.NotNull;

/**
*
* @author Nicholas Roehner
*/
@NoArgsConstructor
public class Sample extends SharableObjBase {

    @Getter
    @Setter
    @Reference
    protected BioDesign bioDesign;

    @NotNull
    @Getter
    @Setter
    @Reference
    protected Container container;
    
    @Getter
    protected Set<Parameter> parameters;

    @Getter
    @Setter
    @ReferenceCollection
    protected Set<Sample> parentSamples;

    public Sample(String name, Container container, Person author) {
        super(name, author);
        this.container = container;
    }

    public Sample(String name, String description, Container container, Person author) {
        super(name, author, description);
        this.container = container;
    }
    
    public Parameter createParameter(double value, Variable variable) {
        Parameter parameter = new Parameter(value, variable);
        addParameter(parameter);
        return parameter;
    }
    
    public void addParameter(Parameter parameter) {
    	if (parameters == null) {
    		parameters = new HashSet<Parameter>();
    	}
    	parameters.add(parameter);
    }
    
    public void addParentSample(Sample parentSample) {
    	if (parentSamples == null) {
    		parentSamples = new HashSet<Sample>();
    	}
    	parentSamples.add(parentSample);
    }

}

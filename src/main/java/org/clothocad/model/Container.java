package org.clothocad.model;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.Reference;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

/**
*
* @author Nicholas Roehner
*/
@NoArgsConstructor
public class Container extends SharableObjBase {

    @Getter
    @Setter
    @Reference
    protected Container container;
    
    @Getter
    @Setter
    protected List<Integer> coordinates;
    
    @Getter
    protected Set<Parameter> parameters;
    
    @Getter
    @Setter
    protected ContainerType type;

    public Container(String name, Person author) {
        super(name, author);
    }

    public Container(String name, String description, Person author) {
        super(name, author, description);
    }
    
    public void addCoordinate(Integer coordinate) {
    	if (coordinates == null) {
    		coordinates = new ArrayList<Integer>();
    	}
    	coordinates.add(coordinate);
    }
    
    public Parameter createParameter(String name, double value, String variable, String units) {
        Parameter parameter = new Parameter(name, value, variable, units);
        addParameter(parameter);
        return parameter;
    }
    
    public void addParameter(Parameter parameter) {
    	if (parameters == null) {
    		parameters = new HashSet<Parameter>();
    	}
    	parameters.add(parameter);
    }
    
    // Feel free to add more of these
    public static enum ContainerType {
        BEAKER, BOX, FLASK, FRIDGE, INCUBATOR, PLATE, RACK, TUBE, WELL;
    }

}

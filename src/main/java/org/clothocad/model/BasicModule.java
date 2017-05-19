package org.clothocad.model;

import org.clothocad.core.persistence.annotations.ReferenceCollection;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import java.util.HashSet;
import java.util.Set;

import javax.validation.constraints.NotNull;
import javax.validation.constraints.Size;

/**
*
* @author Nicholas Roehner
*/
@NoArgsConstructor
public class BasicModule extends Module {

    @NotNull
    @Size(min=1)
    @Getter
    @Setter
    @ReferenceCollection
    protected Set<Feature> features;

    public BasicModule(String name, String role, Set<Feature> features, Person author) {
        super(name, role, author);
        this.features = features;
    }

    public BasicModule(String name, String description, String role, Set<Feature> features, Person author) {
        super(name, description, role, author);
        this.features = features;
    }

    public BasicModule(String name, String role, Person auth) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    @SuppressWarnings("Convert2Diamond")
    public void addFeature(Feature feature) {
    	if (features == null) {
    		features = new HashSet<Feature>();
    	}
    	features.add(feature);
    }

}

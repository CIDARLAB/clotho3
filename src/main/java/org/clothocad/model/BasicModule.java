package org.clothocad.model;

import java.util.Set;

import javax.validation.constraints.NotNull;

import org.clothocad.core.persistence.annotations.ReferenceCollection;

import lombok.Getter;
import lombok.Setter;

public class BasicModule extends Module {
	
	@NotNull
	@Getter
	@Setter
	@ReferenceCollection
	private Set<Feature> features;
	
	public BasicModule(String name, ModuleRole role, Set<Feature> features, Person author) {
		super(name, role, author);
		this.features = features;
	}
	
	public BasicModule(String name, ModuleRole role, Set<Feature> features, Set<Influence> influences, Person author) {
		super(name, role, influences, author);
		this.features = features;
	}
	
	public BasicModule(String name, String description, ModuleRole role, Set<Feature> features, Person author) {
		super(name, description, role, author);
		this.features = features;
	}
	
	public BasicModule(String name, String description, ModuleRole role, Set<Feature> features, Set<Influence> influences, Person author) {
		super(name, description, role, influences, author);
		this.features = features;
	}
	
	public void addSubModule(Feature feature) {
		features.add(feature);
	}
	
	public void removeSubModule(Feature feature) {
		features.remove(feature);
	}
	
}

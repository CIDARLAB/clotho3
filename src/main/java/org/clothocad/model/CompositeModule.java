package org.clothocad.model;

import java.util.Set;

import javax.validation.constraints.NotNull;

import org.clothocad.core.persistence.annotations.ReferenceCollection;

import lombok.Getter;
import lombok.Setter;

public class CompositeModule extends Module {

	@NotNull
	@Getter
	@Setter
	@ReferenceCollection
	private Set<Module> subModules;
	
	public CompositeModule(String name, ModuleRole role, Set<Module> subModules, Person author) {
		super(name, role, author);
		this.subModules = subModules;
	}
	
	public CompositeModule(String name, ModuleRole role, Set<Module> subModules, Set<Influence> influences, Person author) {
		super(name, role, influences, author);
		this.subModules = subModules;
	}
	
	public CompositeModule(String name, String description, ModuleRole role, Set<Module> subModules, Person author) {
		super(name, description, role, author);
		this.subModules = subModules;
	}
	
	public CompositeModule(String name, String description, ModuleRole role, Set<Module> subModules, Set<Influence> influences, Person author) {
		super(name, description, role, influences, author);
		this.subModules = subModules;
	}
	
	public void addSubModule(Module subModule) {
		subModules.add(subModule);
	}
	
	public void removeSubModule(Module subModule) {
		subModules.remove(subModule);
	}
	
}

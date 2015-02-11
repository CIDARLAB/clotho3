package org.clothocad.model;

import java.util.Set;

import javax.validation.constraints.NotNull;
import javax.validation.constraints.Size;

import org.clothocad.core.persistence.annotations.ReferenceCollection;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

/**
*
* @author Nicholas Roehner
*/
@NoArgsConstructor
public class CompositeModule extends Module {

	@NotNull
	@Size(min=2)
	@Getter
	@Setter
	@ReferenceCollection
	private Set<Module> subModules;
	
	public CompositeModule(String name, ModuleRole role, Set<Module> subModules, Person author) {
		super(name, role, author);
		this.subModules = subModules;
	}
	
	public CompositeModule(String name, String description, ModuleRole role, Set<Module> subModules, Person author) {
		super(name, description, role, author);
		this.subModules = subModules;
	}
	
}

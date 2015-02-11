package org.clothocad.model;

import java.util.LinkedList;
import java.util.List;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import org.clothocad.core.persistence.annotations.ReferenceCollection;

/**
*
* @author Nicholas Roehner
*/
@NoArgsConstructor
public class Assembly {
	
	@Getter
	@Setter
	@ReferenceCollection
	private List<Part> parts = new LinkedList<Part>();
	
	@Getter
	@Setter
	private List<Assembly> subAssemblies = new LinkedList<Assembly>();
	
	public Assembly(Part parentPart) {
		parentPart.getAssemblies().add(this);
	}
	
	public Assembly(Assembly parentAssembly) {
		parentAssembly.getSubAssemblies().add(this);
	}

}

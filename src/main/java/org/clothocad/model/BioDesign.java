package org.clothocad.model;

import java.util.HashSet;
import java.util.Set;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.Reference;
import org.clothocad.core.persistence.annotations.ReferenceCollection;

/**
*
* @author Nicholas Roehner
*/
@NoArgsConstructor
public class BioDesign extends SharableObjBase {
	
	@Getter
	@Setter
	@Reference
	private Module module;
	
	@Getter
	@Setter
	private Set<Parameter> parameters = new HashSet<Parameter>();
	
	@Getter
	@Setter
	@ReferenceCollection
	private Set<Part> parts = new HashSet<Part>();
	
	@Getter
	@Setter
	@ReferenceCollection
	private Set<Polynucleotide> polynucleotides = new HashSet<Polynucleotide>();
	
	@Getter
	@Setter
	@ReferenceCollection
	private Set<Strain> strains = new HashSet<Strain>();
	
	@Getter
	@Setter
	@ReferenceCollection
	private Set<Medium> media = new HashSet<Medium>();
	
	@Getter
	@Setter
	@ReferenceCollection
	private Set<BioDesign> subDesigns = new HashSet<BioDesign>();
	
	@Getter
	@Setter
	@Reference
	private BioDesign parentDesign;
	
	public BioDesign(String name, String description, Person author) {
		super(name, author, description);
	}
	
	public BioDesign(String name, Person author) {
		super(name, author);
	}
		
}

package org.clothocad.model;

import java.util.HashSet;
import java.util.Set;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.ReferenceCollection;

/**
 * 
 * @author Nicholas Roehner
 */ 
@NoArgsConstructor
public class Strain extends SharableObjBase {
	
	@Getter
	@Setter
	@ReferenceCollection
	private Set<Polynucleotide> polynucleotides = new HashSet<Polynucleotide>();
	
	@Getter
	@Setter
	private Strain parentStrain;
	
	public Strain(String name, Person author) {
		super(name, author);
	}
	
	public Strain(String name, String description, Person author) {
		super(name, author, description);
	}

}

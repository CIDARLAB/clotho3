package org.clothocad.model;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import org.clothocad.core.datums.SharableObjBase;

/**
 * 
 * @author Nicholas Roehner
 */ 
@NoArgsConstructor
public class Medium extends SharableObjBase {
	
	@Getter
	@Setter
	private Strain parentMedium;
	
	public Medium(String name, Person author) {
		super(name, author);
	}
	
	public Medium(String name, String description, Person author) {
		super(name, author, description);
	}

}

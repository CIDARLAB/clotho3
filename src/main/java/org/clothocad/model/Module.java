package org.clothocad.model;

import java.util.HashSet;
import java.util.Set;

import javax.validation.constraints.NotNull;

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
public abstract class Module extends SharableObjBase {
	
	@NotNull
	@Getter
	@Setter
	protected ModuleRole role;
	
	@Getter
	@Setter
	@ReferenceCollection
	protected Set<Influence> influences = new HashSet<Influence>();
	
	@Getter
	@Setter
	protected Module parentModule;
	
	protected Module(String name, ModuleRole role, Person author) {
		super(name, author);
		this.role = role;
	}
	
	protected Module(String name, String description, ModuleRole role, Person author) {
		super(name, author, description);
		this.role = role;
	}
	
	// Feel free to add more of these
    public static enum ModuleRole {
    	TRANSCRIPTION, TRANSLATION, EXPRESSION;
    }

}

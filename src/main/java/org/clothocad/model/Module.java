package org.clothocad.model;

import java.util.HashSet;
import java.util.Set;

import javax.validation.constraints.NotNull;

import lombok.Getter;
import lombok.Setter;

import org.clothocad.core.datums.SharableObjBase;
import org.clothocad.core.persistence.annotations.ReferenceCollection;

public class Module extends SharableObjBase {
	
	@NotNull
	@Getter
	@Setter
	protected ModuleRole role;
	
	@Getter
	@Setter
	@ReferenceCollection
	protected Set<Influence> influences;
	
	protected Module(String name, ModuleRole role, Person author) {
		super(name, author);
		this.role = role;
		influences = new HashSet<Influence>();
	}
	
	protected Module(String name, ModuleRole role, Set<Influence> influences, Person author) {
		super(name, author);
		this.role = role;
		this.influences = influences;
	}
	
	protected Module(String name, String description, ModuleRole role, Person author) {
		super(name, author, description);
		this.role = role;
		influences = new HashSet<Influence>();
	}
	
	protected Module(String name, String description, ModuleRole role, Set<Influence> influences, Person author) {
		super(name, author, description);
		this.role = role;
		this.influences = influences;
	}
	
	public void addInfluence(Influence influence) {
		influences.add(influence);
	}
	
	public void removeInfluence(Influence influence) {
		influences.remove(influence);
	}
	
	// Feel free to add more of these
    public static enum ModuleRole {
    	TRANSCRIPTION, TRANSLATION, EXPRESSION;
    }

}

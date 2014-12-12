package org.clothocad.model;

import javax.validation.constraints.NotNull;

import lombok.Getter;
import lombok.Setter;

import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.persistence.annotations.Reference;

public class Parameter extends ObjBase {

	@NotNull
	@Getter
	@Setter
	private double value;
	
	@NotNull
	@Getter
	@Setter
	@Reference
	private Variable variable;
	
	@NotNull
	@Getter
	@Setter
	@Reference
	private Units units;
	
	public Parameter(String name, double value, Variable variable, Units units) {
		super(name);
		this.value = value;
		this.variable = variable;
		this.units = units;
	}
}

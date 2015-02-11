package org.clothocad.model;

import javax.validation.constraints.NotNull;

import lombok.Getter;
import lombok.Setter;

import org.clothocad.core.persistence.annotations.Reference;

public class Parameter {

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
	
	public Parameter(BioDesign design, double value, Variable variable, Units units) {
		design.getParameters().add(this);
		this.value = value;
		this.variable = variable;
		this.units = units;
	}
}

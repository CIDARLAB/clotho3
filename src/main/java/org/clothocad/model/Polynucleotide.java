package org.clothocad.model;

import java.io.Serializable;
import java.util.List;
import java.util.Date;

import javax.validation.constraints.NotNull;

import lombok.Getter;
import lombok.Setter;
import lombok.NoArgsConstructor;

import org.clothocad.core.datums.SharableObjBase;

@NoArgsConstructor
public class Polynucleotide extends SharableObjBase implements Serializable {

	@Getter
	@Setter
	private String accession;
	
	@NotNull
	@Getter
	@Setter
	private boolean isLinear, isSingleStranded;
	
	@Getter
	@Setter
	private Date submissionDate;

	@Getter
	@Setter
	private List<Highlight> highlights;

	@Getter
	@Setter
	private String sequence;
	
	@Getter
	@Setter
	private Polynucleotide parentPolynucleotide;
	
	public Polynucleotide(String name, boolean isLinear, boolean isSingleStranded, Person author) {
		super(name, author);
		this.isLinear = isLinear;
		this.isSingleStranded = isSingleStranded;
	}
	
	public Polynucleotide(String name, String description, boolean isLinear, boolean isSingleStranded, 
			Person author) {
		super(name, author, description);
		this.isLinear = isLinear;
		this.isSingleStranded = isSingleStranded;
	}

}

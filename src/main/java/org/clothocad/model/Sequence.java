package org.clothocad.model;

import java.util.HashSet;
import java.util.Set;

import javax.validation.constraints.NotNull;
import javax.validation.constraints.Pattern;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import org.clothocad.core.datums.ObjBase;

/**
* @author nroehner
*/
@NoArgsConstructor
public abstract class Sequence extends ObjBase {

	@NotNull
    @Getter
    @Setter
    @Pattern(regexp="[ATUCGRYKMSWBDHVN]*", flags={Pattern.Flag.CASE_INSENSITIVE})
	protected String sequence;
	
	@Getter
	@Setter
	protected Set<Annotation> annotations;
	
	protected Sequence(String name, String sequence) {
		super(name);
		this.sequence = sequence;
		annotations = new HashSet<Annotation>();
	}
	
	protected Sequence(String name, String sequence, Set<Annotation> annotations) {
		super(name);
		this.sequence = sequence;
		this.annotations = annotations;
	}
	
	public void addAnnotation(Annotation anno) {
		annotations.add(anno);
	}
	
}

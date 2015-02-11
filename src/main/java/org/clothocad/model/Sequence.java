package org.clothocad.model;

import java.util.LinkedList;
import java.util.List;

import javax.validation.constraints.NotNull;
import javax.validation.constraints.Pattern;

import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

import org.clothocad.core.datums.SharableObjBase;

/**
* @author Nicholas Roehner
*/
@NoArgsConstructor
public abstract class Sequence extends SharableObjBase {

	@NotNull
    @Getter
    @Setter
    @Pattern(regexp="[ATUCGRYKMSWBDHVN]*", flags={Pattern.Flag.CASE_INSENSITIVE})
	protected String sequence;
	
	@Getter
	@Setter
	protected List<Annotation> annotations = new LinkedList<Annotation>();
	
	protected Sequence(String name, String sequence, Person author) {
		super(name, author);
		this.sequence = sequence;
	}
	
	protected Sequence(String name, String description, String sequence, Person author) {
		super(name, author, description);
		this.sequence = sequence;
	}
	
}

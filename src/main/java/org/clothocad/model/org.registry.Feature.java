package org.clothocad.model;

import java.util.List;

import lombok.Data;
import lombok.NoArgsConstructor;
//import org.clothocad.core.persistence.annotations.Reference;

@Data()
@NoArgsConstructor
public class org.registry.Feature {

	private boolean direction;
	private String title, type, id;
	private int startpos, endpos;

}
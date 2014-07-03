package org.clothocad.model;

import lombok.Data;
import lombok.NoArgsConstructor;

@Data()
@NoArgsConstructor
public class org.registry.Feature {

	private boolean direction;
	private String title, type, id;
	private int startpos, endpos;

}
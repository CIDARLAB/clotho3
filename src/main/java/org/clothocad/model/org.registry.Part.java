package org.clothocad.model;

import java.util.List;

import lombok.Data;
import lombok.NoArgsConstructor;

@Data()
@NoArgsConstructor
public class org.registry.Part {

	private List<Feature> features = new List<Feature>();
	private List<Subscar> specified_subscars = new List<Subscar>();
	private List<Subpart> specified_subparts = new List<Subpart>();
	private List<Subpart> deep_subparts = new List<Subpart>();
	private List<String> sequences = List<String>();
	private List<String> twins = List<String>();
	

	private String release_status, references, part_nickname, parameters,
		part_url, part_type, sample_status, part_results, samples,
		part_short_name, part_rating, part_id, part_short_desc, categories, groups,
		part_name, part_author, part_entered;
}
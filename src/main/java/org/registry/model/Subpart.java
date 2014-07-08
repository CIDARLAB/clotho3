package org.registry.model;



import lombok.Data;
import lombok.NoArgsConstructor;
import org.clothocad.core.datums.ObjBase;

@Data()
@NoArgsConstructor
public class Subpart extends ObjBase implements Subscar {

	private String part_type, part_nickname, part_name, part_id, part_short_desc;

}
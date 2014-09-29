package org.registry.model;

import java.util.ArrayList;
import java.util.List;

import lombok.Data;
import lombok.NoArgsConstructor;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.SharableObjBase;

@Data()
@NoArgsConstructor
public class Enzyme extends SharableObjBase {

	private List<Polypeptide> polypeptides = new ArrayList<Polypeptide>();

	private String id, version, name, description;

}
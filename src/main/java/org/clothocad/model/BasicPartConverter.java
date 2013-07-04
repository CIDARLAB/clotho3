/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.model;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.schema.Converter;
import org.clothocad.core.schema.InferredSchema;
import org.clothocad.core.schema.Schema;

/**
 *
 * @author spaige
 */
public class BasicPartConverter extends Converter<BasicPart> {
        static Set<String> names = new HashSet<String>();
        static {
            names.add("eugene.dom.components.Part");
        }
    public BasicPartConverter(Persistor p) {
        super(p.get(Schema.class, p.resolveSelector("BasicPart", false)), new HashSet<Schema>(), names);
    }

    @Override
    protected BasicPart guardedConvert(Map data, String schemaName) {
        switch (schemaName){
            case "eugene.dom.components.Part":
                return convertEugenePartToBasicPart(data);
            default:
                return null;
        }
        
    }

    public static BasicPart convertEugenePartToBasicPart(Map<String, Object> eugenePart) {
        BasicPart part = new BasicPart(eugenePart.get("name").toString(), null, eugenePart.get("sequence").toString(), new FreeForm(), null);

        try {
            Part.PartFunction function = Part.PartFunction.valueOf(eugenePart.get("partType").toString().toUpperCase());
            part.setType(function);
        } catch (IllegalArgumentException e) {
        }

        return part;
    }

    @Override
    protected BasicPart guardedConvert(Map data, Schema type) {
        if (type instanceof InferredSchema ){
            return guardedConvert(data, type.getName());
        }
        
        throw new UnsupportedOperationException("No class schemas supported"); //To change body of generated methods, choose Tools | Templates.
    }
}

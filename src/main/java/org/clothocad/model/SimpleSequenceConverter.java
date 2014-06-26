/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.clothocad.model;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.schema.Converter;
import org.clothocad.core.schema.InferredSchema;
import org.clothocad.core.schema.Schema;
import static org.clothocad.model.NucSeqConverter.convertSimpleSeqToNucSeq;
import static org.clothocad.model.NucSeqConverter.names;

/**
 *
 * @author prashantvaidyanathan
 */
public class SimpleSequenceConverter extends Converter<SimpleSequence>
{

    static Set<String> names = new HashSet<String>(); // List of Schemas that can be converted into SimpleSequence
    static  
    {
        names.add("org.clothocad.model.NucSeq");
    }
    
    public SimpleSequenceConverter(Persistor p) 
    {
        super(p.get(Schema.class, new ObjectId("org.clothocad.model.SimpleSequence")), new HashSet<Schema>(), names);
    }
    
    @Override
    protected SimpleSequence guardedConvert(Map data, String schemaName) 
    {
        switch(schemaName)
        {
            case "org.clothocad.model.NucSeq":      //If Schema is of Type Simple Seq
                return convertNucSeqToSimpleSeq(data);
            default:
                return null;
        }
    }
    
    public static SimpleSequence convertNucSeqToSimpleSeq(Map<String,Object> nSeq)
    {
        
        SimpleSequence simpleseq = new SimpleSequence(nSeq.get("name").toString(),nSeq.get("sequence").toString()); //Invoke the NucSeq Constructor that creates an object with the Sequence as the input argument. 
        if(nSeq.containsKey("_id"))
        {
            simpleseq.setId(new ObjectId(nSeq.get("_id").toString()));
        }
        return simpleseq;
    }
    
    
    @Override
    protected SimpleSequence guardedConvert(Map data, Schema type) 
    {
        if(type instanceof InferredSchema)
        {
            return guardedConvert(data,type.getName());
        }
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}

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

/**
 *
 * @author prashantvaidyanathan
 */
public class NucSeqConverter extends Converter<NucSeq>
{
    
    static Set<String> names = new HashSet<String>(); // List of Schemas that can be converted into NucSeq
    static  
    {
        names.add("org.clothocad.model.SimpleSequence");
    }
    public NucSeqConverter(Persistor p) 
    {
        super(p.get(Schema.class, new ObjectId("org.clothocad.model.NucSeq")), new HashSet<Schema>(), names);
    }

    
    @Override
    protected NucSeq guardedConvert(Map data, String schemaName) 
    {
        switch(schemaName)
        {
            case "org.clothocad.model.SimpleSeqeuence":      //If Schema is of Type Simple Seq
                return convertSimpleSeqToNucSeq(data);
            default:
                return null;
        }
    }
    
    public static NucSeq convertSimpleSeqToNucSeq(Map<String,Object> simpleSeq)
    {
        NucSeq nseq = new NucSeq(simpleSeq.get("sequence").toString()); //Invoke the NucSeq Constructor that creates an object with the Sequence as the input argument. 
        if(simpleSeq.containsKey("_id"))
        {
            nseq.setId(new ObjectId(simpleSeq.get("_id").toString()));
        }
        return nseq;
    }
    
    @Override
    protected NucSeq guardedConvert(Map data, Schema type) {
        
        if(type instanceof InferredSchema)
        {
            return guardedConvert(data,type.getName());
        }
        
        
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}

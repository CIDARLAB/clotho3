/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.clothocad.model;

import org.clothocad.core.execution.SequenceConverters;
import java.lang.reflect.Method;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.schema.Converter;
import org.clothocad.core.schema.InferredSchema;
import org.clothocad.core.schema.RunScripts;
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
            case "org.clothocad.model.SimpleSequence":      //If Schema is of Type Simple Seq
                try 
                {
                    return convertSimpleSeqToNucSeq(data);
                } 
                catch (Exception ex) 
                {
                    Logger.getLogger(NucSeqConverter.class.getName()).log(Level.SEVERE, null, ex);
                }
            default:
                return null;
        }
    }
    
    public static NucSeq convertSimpleSeqToNucSeq(Map<String,Object> simpleSeq) throws Exception
    {
        
        String ConvertedSequence = "";
        
        Map<Object,Method> functionMap = new HashMap<Object,Method>();
        HashMap<Object,Map<Object,Method>> conversionMap = new HashMap<Object,Map<Object,Method>>();
        
        
        //String inpSequence = 
        Class[] parameterTypes = new Class[1];
        //parameterTypes[0] = simpleSeq.get("sequence").toString().getClass();
        parameterTypes[0] = simpleSeq.get("sequence").getClass();
        Method getSequence = SequenceConverters.class.getMethod("getSequence", parameterTypes);
        
        functionMap.put("sequence", getSequence);
        conversionMap.put("sequence", functionMap);
        
        RunScripts mapper = new RunScripts();
        String res ="";
        res = (String)mapper.runningfunction(mapper, conversionMap.get("sequence").get("sequence"), simpleSeq.get("sequence").toString());
        
        System.out.println("This is the result of the function :"+res);
        
        NucSeq nseq = new NucSeq(res, (Person) simpleSeq.get("author")); 
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

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.clothocad.core.schema;

import java.lang.reflect.Field;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.model.NucSeq;
import org.clothocad.model.SimpleSequence;

/**
 *
 * @author prashantvaidyanathan
 */
public class SequenceConverters extends Converter{

    
    public static Persistor p;
    public SimpleSequence naiveConvertNucSeqToSimpleSeq(NucSeq nucObj)
    {
        SimpleSequence simpleSeq = new SimpleSequence(nucObj.getName(),nucObj.getSeq());
        return simpleSeq;
    }
    
    public NucSeq naiveConvertNucSeqToSimpleSeq(SimpleSequence simpleObj)
    {
        NucSeq nucSeq = new NucSeq(simpleObj.getSequence());
        return nucSeq;
    }
    //resolveSchemaFromClassName
    
    public static SimpleSequence convertNucSeqToSimpleSeq(NucSeq nucObj) throws IllegalArgumentException, IllegalAccessException
    {
        
        Schema nucSchema = new InferredSchema("nucseq");
        Schema simpleSchema = new InferredSchema("simpleseq");
        Map <String, Object> nucMap = new HashMap<String, Object>();
        
        Field fld[] = NucSeq.class.getDeclaredFields();
        for(Field  xfield: fld)
        {
            xfield.setAccessible(true);
            nucMap.put(xfield.getName(), xfield.get(nucObj));
            //System.out.println(xfield.getName());
            //System.out.println(xfield.get(nucObj).toString());
        }
        
        
        //Need something like Converter<SimpleSeq>
        //SimpleSequence simpleseqObj = converter.convert(nucMap, nucSchema);
        
        
        
        return null;
    
    }
    
    
    public SequenceConverters(Schema convertsTo, Set canConvert, Set canConvertNames) {
        super(convertsTo, canConvert, canConvertNames);
    }

    @Override
    protected ObjBase guardedConvert(Map data, Schema type) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    protected ObjBase guardedConvert(Map data, String schemaName) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
}

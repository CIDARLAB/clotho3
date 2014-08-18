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
import javafx.util.Pair;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.model.NucSeq;
import org.clothocad.model.SimpleSequence;

/**
 *
 * @author prashantvaidyanathan
 */
public class ConverterMapper 
{
    /*
    Schema nucSchema = new InferredSchema("nucseq");
        Schema simpleSchema = new InferredSchema("simpleseq");
        Map <String, Object> nucMap = new HashMap<String, Object>();
        
        Field fld[] = NucSeq.class.getDeclaredFields();
        for(Field  xfield: fld)
        {
            xfield.setAccessible(true);
            nucMap.put(xfield.getName(), xfield.get(nucObj));
        }
    */
    
    public static void testclass(Object seqobj,String propname) throws Exception
    {
        Field xfield;
        seqobj.getClass().getDeclaredField(propname).setAccessible(true);
        xfield = seqobj.getClass().getDeclaredField(propname);
        xfield.setAccessible(true);
        System.out.println("Property "+propname + ", of Object "+ seqobj.getClass() +  " : "+ xfield.get(seqobj));
    }
    
}

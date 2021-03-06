/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.communication;

import com.google.common.collect.Lists;
import java.util.ArrayList;
import java.util.List;
import org.clothocad.core.datums.Argument;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.datums.ObjectId;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author spaige
 */
public class RunTest extends AbstractServerAPITest  {


    @Test 
    public void resolveIdToObject(){
        Argument[] functionArguments = {new Argument("testArgument", ObjBase.class)};
        List<Object> resolvedIds = api.resolveIds(functionArguments, Lists.asList(ids.get(0).toString(), new Object[0]));
        assertEquals(resolvedIds.get(0), persistor.getAsJSON(ids.get(0)));
    }
    
    @Test
    public void resolveNonExistentIdToObject(){
        Argument[] functionArguments = {new Argument("testArgument", ObjBase.class)};
        List<Object> resolvedIds = api.resolveIds(functionArguments, Lists.asList(new ObjectId().toString(), new Object[0])); 
        assertTrue(resolvedIds.get(0) instanceof String);
    }
    
    @Test
    public void doNotResolveIdToObject(){
        Argument[] functionArguments = {new Argument("testArgument", ObjectId.class)};
        List<Object> resolvedIds = api.resolveIds(functionArguments, Lists.asList(ids.get(0).toString(), new Object[0]));
        assertTrue(resolvedIds.get(0) instanceof String);
    }
    
    @Test 
    public void convertSomeToObject(){
        Argument[] functionArguments = {new Argument("testArgument", ObjectId.class), 
                                        new Argument("anotherTestArgument", ObjBase.class),
                                        new Argument("thirdTestArgument", Number.class)};
        List<Object> arguments = new ArrayList();
        arguments.add(ids.get(0).toString());
        arguments.add(ids.get(0).toString());
        arguments.add(3);
        List<Object> resolvedIds = api.resolveIds(functionArguments, arguments);
        assertTrue(resolvedIds.get(0) instanceof String);
        assertEquals(resolvedIds.get(1), persistor.getAsJSON(ids.get(0)));
        assertEquals(resolvedIds.get(2), 3);
        
    }
}

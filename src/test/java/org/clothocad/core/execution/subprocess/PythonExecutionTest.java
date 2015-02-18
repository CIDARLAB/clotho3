/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.execution.subprocess;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.clothocad.core.communication.AbstractServerAPITest;
import org.junit.Test;
import static org.junit.Assert.*;
/**
 *
 * @author spaige
 */
public class PythonExecutionTest extends AbstractServerAPITest {
    @Test
    public void canAccessAttrs(){
        final Map<String, Object> argument = new HashMap<>();
        argument.put("accessMe", "passed!");
        final List<Object> args = new ArrayList<>();
        args.add(argument);
        assertEquals(pythonSubprocessExec("def run(obj): return obj.accessMe", args), "passed!");
    }
    
    @Test
    public void canSetNewAttrs(){
        final Map<String,Object> argument = new HashMap<>();
        final List<Object> args = new ArrayList<>();
        args.add(argument);
        assertEquals(pythonSubprocessExec("def run(obj):\n"
                + "    obj.newAttr = 5\n"
                + "    return obj.newAttr\n", args), 5);
    }
    
    @Test
    public void testPythonExecution(){
        assertEquals(pythonSubprocessExec("def run(): return 'passed!'", new ArrayList<>()), "passed!");
    }
   
    @Test
    public void testReturnList(){
        List result = (List) pythonSubprocessExec(
                "def run(): return [1,2,3]", new ArrayList<>());
        assertEquals(result, Arrays.asList(1,2,3));
    }
    
    @Test
    public void testReturnClass(){
        Map<String, Object> result = (Map) pythonSubprocessExec(
                  "class TestClass(object):\n"
                + "    def __init__(self):\n"
                + "        self.test = 'passed!'\n"
                + "\n"
                + "def run(): return TestClass()", new ArrayList<>());
        
        assertEquals(result.get("test"), "passed!");
    }

    @Test
    public void testGetAPI(){
        Map<String,Object> result = (Map) pythonSubprocessExec(
                "def run(): return clotho.get('org.clothocad.core.datums.Module')",
                new ArrayList<>());
        
        assertEquals(result.get("name"), "Module");
    }
    
    private Object pythonSubprocessExec(String code, List<Object> args){
        final Map<String, Object> sourceJSON = new HashMap<>();
        sourceJSON.put("code", code);
        sourceJSON.put("name", "tester");
        sourceJSON.put("language", "python");
        return SubprocessExec.run(api, sourceJSON, args, new SubprocessExec.EventHandler() {

            @Override
            public void onFail(byte[] standardError) {
                throw new RuntimeException(new String(standardError, StandardCharsets.UTF_8));
            }

            @Override
            public void onSuccess(byte[] standardError) {
            }
                    
                });
    }
    
   /* public void subprocessExecutePCRPrediction(){
        final Map<String,Object> function = loadJSON(Paths.get("src", "test", "resources", "authoredJSON", "clotho.functions.dna.pcr.json"));
        List<Object> args = function.get(args);
    } 
    
    
    private static Map<String,Object>loadJSON(Path file){
        
    }*/
}

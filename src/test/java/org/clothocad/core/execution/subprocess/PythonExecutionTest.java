/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.execution.subprocess;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.junit.Test;
import static org.junit.Assert.*;
/**
 *
 * @author spaige
 */
public class PythonExecutionTest {
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
    
    private Object pythonSubprocessExec(String code, List<Object> args){
        final Map<String, Object> sourceJSON = new HashMap<>();
        sourceJSON.put("code", code);
        sourceJSON.put("name", "tester");
        sourceJSON.put("language", "python");
        return SubprocessExec.run(null, sourceJSON, args, new SubprocessExec.EventHandler() {

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

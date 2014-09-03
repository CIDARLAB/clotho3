/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.execution;

import com.google.inject.Injector;
import javax.script.ScriptException;
import org.clothocad.core.communication.ServerSideAPI;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.security.ClothoRealm;
import org.clothocad.core.util.TestUtils;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author spaige
 */
public class ScriptEngineTest {

    MetaEngine engine;
    
    public ScriptAPI genAPI(){
        Injector injector = TestUtils.getDefaultTestInjector();
        return new ScriptAPI(new ServerSideAPI(null, injector.getInstance(Persistor.class), null, injector.getInstance(ClothoRealm.class), null));
    }
    
    @Before
    public void setup(){
        engine = new MetaEngine();
    }
    
    @Test
    public void testJavaScriptAPI() throws ScriptException{
        Object result = engine.eval("clotho.newProperty = 5; clotho.newProperty;", Language.JAVASCRIPT, genAPI());
        assertEquals(5, result);
        result = engine.eval("clotho.newProperty", Language.JAVASCRIPT, genAPI());
        assertEquals(5, result);
    }
}

package org.clothocad.core.testers;

import org.junit.Before;
import org.junit.Test;
import org.clothocad.core.layers.communication.mind.Mind;
import org.clothocad.core.util.FileUtils;
import org.junit.After;

public class MindTest {
    Mind mind = new Mind();
    
    @Before
    public void setUp() {
        System.out.println("setup");
    }
    
    @After
    public void tearDown() {
        System.out.println("tear");
    }

    @Test
    public void testJSScript() throws Exception{
        String script = FileUtils.readFile("mind_test.js");
        assert(mind.runCommand(script));
    }
}

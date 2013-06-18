package org.clothocad.core.testers;



import com.github.jmkgreen.morphia.logging.MorphiaLoggerFactory;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertSame;

import java.io.IOException;
import java.util.Random;

import org.bson.types.ObjectId;
import org.junit.Before;
import org.junit.Test;

import com.mongodb.BasicDBObject;
import java.util.HashMap;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.clothocad.core.aspects.Collector;
import org.clothocad.core.aspects.Persistor;
import org.clothocad.core.datums.ObjBase;
import org.clothocad.core.layers.communication.Router;
import org.clothocad.core.layers.communication.mind.Mind;
import org.clothocad.core.util.FileUtils;
import org.clothocad.model.Institution;
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

package org.clothocad.core.testers;

import java.util.ArrayList;
import java.util.List;
import org.bson.types.ObjectId;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.aspects.Proctor.Module;
import org.clothocad.core.aspects.Proctor.Paver;
import org.clothocad.core.aspects.Proctor.TemplatePaver;
import org.junit.Before;
import org.junit.Test;
import org.clothocad.core.layers.communication.mind.Mind;
import org.clothocad.core.utils.TestUtils;
import org.clothocad.model.Trail;
import org.junit.After;

public class TrailTest {
    Mind mind = new Mind();
    Persistor persistor = TestUtils.getA(Persistor.class);
    
    @Before
    public void setUp() {
        System.out.println("setup");
    }
    
    @After
    public void tearDown() {
        System.out.println("tear");
    }

    @Test
    public void testCreateTrail() throws Exception{
        //construct the contents
        List<Module> contents = new ArrayList<>();
        List<Paver> pavers = new ArrayList<>();
        TemplatePaver tp1 = new TemplatePaver("Introduction", "/app/partials/trail_123_1.html");
        pavers.add(tp1);
        Module module = new Module("The Basics", pavers);
        contents.add(module);
        
        //Make the trail and confirm it is populated
        Trail trail = new Trail("First Biosafety Module", "This is a trail test", contents);
        assert(trail.getContents().get(0).getPavers().get(0).getType().equals("template"));
        
        
        //Save then re-retrieve the trail
        persistor.save(trail);
        ObjectId uuid = trail.getUUID();
        Trail result = persistor.get(Trail.class, uuid);
        assert(result.getTitle().equals("First Biosafety Module"));
        
        
    }
}

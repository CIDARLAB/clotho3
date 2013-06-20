package org.clothocad.core.testers;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.bson.types.ObjectId;
import org.clothocad.core.aspects.Persistor;
import org.clothocad.core.aspects.Proctor.Module;
import org.clothocad.core.aspects.Proctor.Paver;
import org.clothocad.core.aspects.Proctor.TemplatePaver;
import org.clothocad.core.datums.ObjBase;
import org.junit.Before;
import org.junit.Test;
import org.clothocad.core.layers.communication.mind.Mind;
import org.clothocad.model.Instance;
import org.clothocad.model.Trail;
import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;
import org.junit.After;

public class TrailTest {
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
    public void testCreateTrail() throws Exception{
        //construct the contents
        List<Module> contents = new ArrayList<Module>();
        List<Paver> pavers = new ArrayList<Paver>();
        TemplatePaver tp1 = new TemplatePaver("Introduction", "/app/partials/trail_123_1.html");
        pavers.add(tp1);
        Module module = new Module("The Basics", pavers);
        contents.add(module);
        
        //Make the trail and confirm it is populated
        Trail trail = new Trail("First Biosafety Module", "This is a trail test", contents);
        assert(trail.getContents().get(0).getPavers().get(0).getType().equals("template"));
        
        
        //Save then re-retrieve the trail
        Persistor.get().save(trail);
        ObjectId uuid = trail.getUUID();
        Trail result = Persistor.get().get(Trail.class, uuid);
        assert(result.getTitle().equals("First Biosafety Module"));
        
        //var trails = clotho.query({"className":"org.clothocad.model.Trail"});
        //var trail = clotho.get('51c20034507659b64be65a3b');
        //clotho.startTrail('51c20034507659b64be65a3b');
        
//        System.out.println("The end result is:");
//        System.out.println(result.toString());
        
    }
    
    @Test
    public void testInstance() throws Exception {
//            JSONObject json = new JSONObject();
//            
//        //Start of funky jsonobject test
//            json.put("donut", "crispy");
////                JSONArray listy = new JSONArray();
////                listy.put("cindy");
////                listy.put("bob");
////                listy.put("mary");
////                listy.put("tom");
////                    JSONObject inner = new JSONObject();
////                    inner.put("age", "5");
////                listy.put(inner);
////            json.put("listy", listy);
//        //end of funky jsonobject test
//                    
//            System.out.println("The json is : " + json.toString());
//                    
//            Instance instance = new Instance(json);
//            //Save then re-retrieve the trail
//            Persistor.get().save(instance);
//            String uuid = instance.getUUID().toString();
//            Instance ires = (Instance) Persistor.get().get(Instance.class, uuid);
//
//            System.out.println("The isnstance result is:");
//            System.out.println(ires.toString());
    }
}

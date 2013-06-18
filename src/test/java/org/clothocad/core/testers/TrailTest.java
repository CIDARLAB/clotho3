package org.clothocad.core.testers;

import org.bson.types.ObjectId;
import org.clothocad.core.aspects.Collector;
import org.clothocad.core.aspects.Persistor;
import org.junit.Before;
import org.junit.Test;
import org.clothocad.core.layers.communication.mind.Mind;
import org.clothocad.model.Trail;
import org.json.JSONArray;
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
        JSONArray json = new JSONArray("[\n" +
"        {\n" +
"            \"module_title\" : \"The Basics\",\n" +
"            \"pavers\" : [\n" +
"                {\n" +
"                    \"paver_title\" : \"Introduction\",\n" +
"                    \"type\" : \"template\",\n" +
"                    \"template\" : \"/app/partials/trail_123_1.html\"\n" +
"                }\n" +
"            ]\n" +
"        },\n" +
"        {\n" +
"            \"module_title\" : \"Demo Quizzes\",\n" +
"            \"pavers\" : [\n" +
"                {\n" +
"                    \"paver_title\" : \"MC Question\",\n" +
"                    \"type\" : \"quiz\",\n" +
"                    \"content\" : {\n" +
"                        \"uuid\" : \"123456\",\n" +
"                        \"title\" : \"Question 1\",\n" +
"                        \"type\" : \"mc\",\n" +
"                        \"question\" : \"What is a plasmid?\",\n" +
"                        \"options\" : [\n" +
"                            \"A bacterial genome\",\n" +
"                            \"A viral genome\",\n" +
"                            \"Like an organelle\",\n" +
"                            \"None of the above\"\n" +
"                        ]\n" +
"                    }\n" +
"                },\n" +
"                {\n" +
"                    \"paver_title\" : \"MultiPic Question\",\n" +
"                    \"type\" : \"quiz\",\n" +
"                    \"content\" : {\n" +
"                        \"uuid\" : \"67890\",\n" +
"                        \"title\" : \"Choose the right image\",\n" +
"                        \"type\" : \"multipic\",\n" +
"                        \"question\" : \"What does Z-form DNA look like?\",\n" +
"                        \"options\" : [\n" +
"                            {\n" +
"                                \"img_url\" : \"http://www.albany.edu/chemistry/cbb/abc_helix.jpg\",\n" +
"                                \"value\" : \"Like this\"\n" +
"                            },\n" +
"                            {\n" +
"                                \"img_url\" : \"https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcTGSGiPIX0EdI25CuzD1rt384QeNLVHZ4KLMtbYdTZMq0fC6Al89w\",\n" +
"                                \"value\" : \"or this\"\n" +
"                            },\n" +
"                            {\n" +
"                                \"img_url\" : \"http://web.chem.ucsb.edu/~molvisual/Img/142A/B_DNA_side.png\",\n" +
"                                \"value\" : \"maybe this\"\n" +
"                            }\n" +
"                        ]\n" +
"                    }\n" +
"                },\n" +
"                {\n" +
"                    \"paver_title\" : \"Fill-in Question\",\n" +
"                    \"type\" : \"quiz\",\n" +
"                    \"content\" : {\n" +
"                        \"uuid\" : \"klasdhf\",\n" +
"                        \"title\" : \"Another Question\",\n" +
"                        \"type\" : \"fillin\",\n" +
"                        \"question\" : \"What is the reverse complement of <code>ACTGCT?</code>\"\n" +
"                    }\n" +
"                },\n" +
"                {\n" +
"                    \"paver_title\" : \"Essay Question\",\n" +
"                    \"type\" : \"quiz\",\n" +
"                    \"content\" : {\n" +
"                        \"uuid\" : \"asdkgjasdg\",\n" +
"                        \"title\" : \"A much longer question\",\n" +
"                        \"type\" : \"essay\",\n" +
"                        \"question\" : \"Describe reverse complementation in great detail.\"\n" +
"                    }\n" +
"                },\n" +
"                {\n" +
"                    \"paver_title\" : \"Drag-Drop Categories\",\n" +
"                    \"type\" : \"quiz\",\n" +
"                    \"content\" : {\n" +
"                        \"uuid\" : \"dragdrop12345\",\n" +
"                        \"title\" : \"Drag Answers to the right section\",\n" +
"                        \"type\" : \"dragdrop\",\n" +
"                        \"question\" : \"What is each type of molecule?\",\n" +
"                        \"options\" : [\n" +
"                            {\n" +
"                                \"type\" : \"text\",\n" +
"                                \"img_url\" : \"http://www.albany.edu/chemistry/cbb/abc_helix.jpg\",\n" +
"                                \"value\" : \"DNA\"\n" +
"                            },\n" +
"                            {\n" +
"                                \"type\" : \"img\",\n" +
"                                \"img_url\" : \"http://web.chem.ucsb.edu/~molvisual/Img/142A/B_DNA_side.png\",\n" +
"                                \"value\" : \"DNA #2\"\n" +
"                            },\n" +
"                            {\n" +
"                                \"type\" : \"text\",\n" +
"                                \"value\" : \"glucose\"\n" +
"                            },\n" +
"                            {\n" +
"                                \"type\" : \"text\",\n" +
"                                \"value\" : \"cellulose\"\n" +
"                            }\n" +
"                        ],\n" +
"                        \"categories\" : [\n" +
"                            \"monomer\",\n" +
"                            \"oligomer\",\n" +
"                            \"polymer\"\n" +
"                        ]\n" +
"                    }\n" +
"                },\n" +
"                {\n" +
"                    \"paver_title\" : \"Multicheck Question\",\n" +
"                    \"type\" : \"quiz\",\n" +
"                    \"content\" : {\n" +
"                        \"uuid\" : \"sfgsufigdjl\",\n" +
"                        \"title\" : \"Select all the correct answers\",\n" +
"                        \"type\" : \"multicheck\",\n" +
"                        \"question\" : \"Which of the following are true of DNA?\",\n" +
"                        \"options\" : [\n" +
"                            \"It contains Uracil\",\n" +
"                            \"High GC content increases melting temperature\",\n" +
"                            \"It is more stable at high temperatures\",\n" +
"                            \"A-form is the most common\"\n" +
"                        ]\n" +
"                    }\n" +
"                },\n" +
"                {\n" +
"                    \"paver_title\" : \"Ranking Question\",\n" +
"                    \"type\" : \"quiz\",\n" +
"                    \"content\" : {\n" +
"                        \"uuid\" : \"s08ergsho\",\n" +
"                        \"title\" : \"Place the options in the correct order\",\n" +
"                        \"type\" : \"ranking\",\n" +
"                        \"question\" : \"Rank the following by molecular weight\",\n" +
"                        \"options\" : [\n" +
"                            \"Uracil\",\n" +
"                            \"Ribose\",\n" +
"                            \"Sucrose\",\n" +
"                            \"Adenosine\",\n" +
"                            \"Thymine\"\n" +
"                        ]\n" +
"                    }\n" +
"                },\n" +
"                {\n" +
"                    \"paver_title\" : \"True/False Question\",\n" +
"                    \"type\" : \"quiz\",\n" +
"                    \"content\" : {\n" +
"                        \"uuid\" : \"sfgsufigdjl\",\n" +
"                        \"title\" : \"Is it true or false?\",\n" +
"                        \"type\" : \"truefalse\",\n" +
"                        \"question\" : \"GC content increases melting temperature\"\n" +
"                    }\n" +
"                },\n" +
"                {\n" +
"                    \"paver_title\" : \"Matching Question\",\n" +
"                    \"type\" : \"quiz\",\n" +
"                    \"content\" : {\n" +
"                        \"uuid\" : \"sasdfas9d\",\n" +
"                        \"title\" : \"Match each\",\n" +
"                        \"type\" : \"match\",\n" +
"                        \"question\" : \"Match each option with its correct match\",\n" +
"                        \"options\" : [\n" +
"                            \"Uracil\",\n" +
"                            \"DNA\",\n" +
"                            \"sucrose\"\n" +
"                        ],\n" +
"                        \"matches\" : [\n" +
"                            \"monomer\",\n" +
"                            \"polymer\",\n" +
"                            \"oligomer\"\n" +
"                        ]\n" +
"                    }\n" +
"                }\n" +
"            ]\n" +
"        }\n" +
"    ]");
        
        Trail trail = new Trail("Trail to knowhere", "This is a trail test", json);
        assert(trail.getTitle().equals("Trail to knowhere"));
        
        JSONArray conts = trail.getContents();
        System.out.println(conts.length());
        JSONObject obj = conts.getJSONObject(0);
        System.out.println("conts.length " + conts.length());
        String title = obj.getString("module_title");
        System.out.println("title " + title);
        assert(title.equals("The Basics"));
        
        //Save then re-retrieve the trail
        Persistor.get().save(trail);
        ObjectId uuid = trail.getUUID();
        Trail result = Persistor.get().get(Trail.class, uuid);
        assert(result.getTitle().equals("Trail to knowhere"));
        
        conts = result.getContents();
        obj = conts.getJSONObject(0);
        title = obj.getString("module_title");
        System.out.println("title " + title);
        assert(title.equals("The Basics"));
    }
}

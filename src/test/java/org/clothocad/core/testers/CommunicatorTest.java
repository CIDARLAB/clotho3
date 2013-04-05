package org.clothocad.core.testers;

import java.util.ArrayList;
import java.util.List;
import org.clothocad.core.aspects.Ambassador.Ambassador;
import org.clothocad.core.aspects.Persistor;
import org.clothocad.core.datums.View;
import org.clothocad.core.datums.objbases.Person;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.datums.util.FieldType;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.datums.util.ServerScript;
import org.clothocad.core.layers.communication.Channel;
import org.clothocad.core.layers.communication.Communicator;
import org.json.JSONObject;

/**
 * Starts web server on http://localhost:1025
 * 
 * var args = []; args[0] = 'static-admin-instance-is-uuid'; clotho.show('CT-sample-view', args, '{}'); 
 *
 * @author Kelvin Li
 */
public class CommunicatorTest {
	
	// In Communicator test we ``push'' a View datum to the client...
	
    public static void main(String[] args) {
    	makeSampleView();
        
        // create a new Communicator instance (following the Singleton pattern)
        //Communicator.get().sendClientMessage("???", Channel.NOTIFICATION, message);
        
        // what should this message contain?
        // - the view object (in JSON)
        // - 
        
       
        // Ambassador.get();

        /***
        ServerSideAPI api = new ServerSideAPI(null, null);
        api.learn("test2", "var args = []; args[0] = 'static-admin-instance-is-uuid'; clotho.show('CT-sample-view', args, '{}'); ");  //IT LEARNS IT, BUT EXECUTION DOESN'T HAPPEN FOR KNOWN REASONS
        
        //Test sending another Clotho a message

        String kevinIp = "10.0.1.5:86";
//        Ambassador.get().sendClothoMessage(kevinIp, "Hi there Kevin, magic word is 'elephant'");
//        System.out.println("completed");
         ***/
    }


    private static void makeSampleView() {
        try {
            ServerScript canUpdate =
                new ServerScript("return true;",
                                 Language.JavaScript);
          
            String html = " <p>\r\n    <label>Email Address\r\n        <input type=\"text\" name=\"email\" id=\"_widget_id_email\" />\r\n      </label>\r\n    </p>\r\n    <p>\r\n      <label>Display Name\r\n        <input type=\"text\" name=\"displayname\" id=\"_widget_id_displayname\" />\r\n      </label>\r\n    </p>\r\n    <p>\r\n      <label>First Name\r\n        <input type=\"text\" name=\"givenname\" id=\"_widget_id_givenname\" />\r\n      </label>\r\n    </p>\r\n    <p>\r\n      <label>Last Name\r\n        <input type=\"text\" name=\"surname\" id=\"_widget_id_surname\" />\r\n      </label>\r\n    </p>\r\n    <p>\r\n      <label>NickName\r\n        <input type=\"text\" name=\"nickname\" id=\"_widget_id_nickname\" />\r\n      </label>\r\n    </p>";
            String onShow = "";
            String onUpdate = "";//alert(JSON.stringify(person));";
            
            List<ClothoField> inputArgs = new ArrayList<ClothoField>();
            inputArgs.add(new ClothoField("person", FieldType.SCHEMA, Person.getSchema().getId(), 3));
            
            View view = View.create(
                         Person.getAdmin(),
                         "Sample View",
                         "A hard-coded view created by CommunicatorTest",
                         inputArgs,
                         canUpdate, 
                         html,
                         onShow,
                         onUpdate);

            //Change the Id
            JSONObject obj = view.toJSON();
            obj.put("id", "CT-sample-view");
            view = View.deserialize(obj.toString());
                        
            view.save();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}




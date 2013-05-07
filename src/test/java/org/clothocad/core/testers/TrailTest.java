package org.clothocad.core.testers;

import java.util.ArrayList;
import java.util.List;
import org.clothocad.core.aspects.Persistor;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.View;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.layers.communication.Communicator;
import org.clothocad.core.schema.Access;
import org.clothocad.model.Person;
import org.json.JSONObject;
import static org.objectweb.asm.Opcodes.*;

/**
 * Starts web server on http://localhost:1025
 * 
 * var args = []; args[0] = 'static-admin-instance-is-uuid'; clotho.show('CT-sample-view', args, '{}'); 
 *
 * @author Kelvin Li
 */
public class TrailTest {
    public static void main(String[] args) {
        makeContentElement();
        makeTrailElement();
        
        Communicator.get();
    }

    private static void makeTrailElement() {
        
    }

    private static void makeContentElement() {
        try {
            //Construct a youtube view
            Function canUpdate = new Function();
                //new Function("", String[]{}, null, null, 
                //    new Script("outputs.put('is_valid', true);", Language.JAVASCRIPT));
          
            String html = " <p>\r\n    <label>Email Address\r\n        <input type=\"text\" name=\"email\" id=\"_widget_id_email\" />\r\n      </label>\r\n    </p>\r\n    <p>\r\n      <label>Display Name\r\n        <input type=\"text\" name=\"displayname\" id=\"_widget_id_displayname\" />\r\n      </label>\r\n    </p>\r\n    <p>\r\n      <label>First Name\r\n        <input type=\"text\" name=\"givenname\" id=\"_widget_id_givenname\" />\r\n      </label>\r\n    </p>\r\n    <p>\r\n      <label>Last Name\r\n        <input type=\"text\" name=\"surname\" id=\"_widget_id_surname\" />\r\n      </label>\r\n    </p>\r\n    <p>\r\n      <label>NickName\r\n        <input type=\"text\" name=\"nickname\" id=\"_widget_id_nickname\" />\r\n      </label>\r\n    </p>";
            String onShow = "";
            String onUpdate = "alert(JSON.stringify(person));";
            
            List<ClothoField> inputArgs = new ArrayList<ClothoField>();
            inputArgs.add(new ClothoField("person", Person.class, null, null, null, true, Access.PUBLIC));
            
            View view = View.create(
                         null,
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
            //view = View.deserialize(obj.toString());
            //XXX: view.save();

        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}




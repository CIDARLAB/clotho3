/*
Copyright (c) 2010 The Regents of the University of California.
All rights reserved.
Permission is hereby granted, without written agreement and without
license or royalty fees, to use, copy, modify, and distribute this
software and its documentation for any purpose, provided that the above
copyright notice and the following two paragraphs appear in all copies
of this software.

IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.

THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
ENHANCEMENTS, OR MODIFICATIONS.
 */

package org.clothocad.core.testers;

import java.util.LinkedHashMap;
import java.util.Map;
import org.clothocad.core.aspects.Persistor;
import org.clothocad.core.datums.Schema;
import org.clothocad.core.datums.View;
import org.clothocad.core.datums.objbases.Person;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.datums.util.FieldType;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.datums.util.ServerScript;
import org.clothocad.core.util.Logger;
import org.json.JSONException;
import org.json.JSONObject;


/**
 * html escape tool: http://www.htmlescape.net/htmlescape_tool.html
 * This test creates the ShowAPerson view
 * @author John Christopher Anderson
 */

public class T5 {
    public static void main(String[] args) {
//        createShowAPersonView();
    }
//    
//    private static View createShowAPersonView() {
//        try {
//            //Start constructing the bits to describe an View
//            ServerScript canshowit = new ServerScript("return true;", Language.JavaScript);
//
//            ServerScript showit = new ServerScript(
//                "&lt;p&gt;I am a widget. I can post things to the log&lt;/p&gt;<br/> &lt;input type=&quot;submit&quot; class=&quot;login_button&quot; name=&quot;example_log&quot; value=&quot;Send to server (ex)&quot; onclick=&quot;say(&amp;quot;This is an example post&amp;quot;)&quot;&gt;",
//                Language.JavaScript
//            );
//
//            ServerScript canupdateit = new ServerScript("", Language.JavaScript);
//            ServerScript updateit = new ServerScript("", Language.JavaScript);
//
//            Schema personSchema = Person.getSchema();
//
//            List<ClothoField> inputArgs = new LinkedHashList<ClothoField>();
//            ClothoField cindy = new ClothoField(FieldType.SCHEMA, personSchema.getId(), 1);
//            inputArgs.put("cindy", cindy);
//
//
//            //Create an View
//            View view = View.create(Person.getAdmin(),
//                    "ShowAPerson",
//                    "Inputs a Person instance and then shows its contents graphically",
//                    canshowit, canupdateit, showit, updateit);
//
//            //Change the Id
//            JSONObject obj = view.toJSON();
//            obj.put("id", "specific-showapersonview-is-uuid");
//            view = View.deserialize(obj.toString());
//            Persistor.get().persistDatum(view);
//            Logger.log(Logger.Level.INFO, view.getId() + "   " + view.getName() + "\n" + view.getDescription() + "\n...was created successfully, all good!");
//            return view;
//        } catch (JSONException ex) {
//            ex.printStackTrace();
//            System.exit(0);
//            return null;
//        }
//    }
}

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

import java.util.ArrayList;
import java.util.List;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.datums.View;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.datums.util.Language;
import org.json.JSONObject;

/**
 * Tests the NameSpace functionality of the Mind
 * Create a Mind, populate it with stuff, clear it, then run something namespaced
 */
public class T7 {

    /**
     * Construct a youtube-viewing widget
     * @param args 
     */
    public static void main(String[] args) {
/*try {
            /*ServerScript canUpdate =
                new ServerScript("outputs.put('is_valid', true);",
                                 Language.JavaScript);
          
            String html = "\t\t\t\t\t<div id=\"video1\">\r\n\t\t\t\t\t\t<div class=\"view\">\r\n\t\t\t\t\t\t\t<iframe width=\"420\" height=\"315\" src=\"http://www.youtube.com/embed/8aTDCZIhcY8\" frameborder=\"0\" allowfullscreen></iframe>\r\n\t\t\t\t\t\t</div>\r\n\t\t\t\t\t\t<div class=\"view padded\">\r\n\t\t\t\t\t\t\t<p>This is content about the video</p>\r\n\t\t\t\t\t\t\t<p>This is more content</p>\r\n\t\t\t\t\t\t</div>\r\n\t\t\t\t\t</div>\r\n\t\t\t\t\t";
            String onShow = "";
            String onUpdate = "";
            
            List<ClothoField> inputArgs = new ArrayList<ClothoField>();
            View view = View.create(
                         Person.getAdmin(),
                         "YouTube viewer",
                         "A lightweight relay of youtube videos",
                         inputArgs,
                         canUpdate, 
                         html,
                         onShow,
                         onUpdate);

            //Change the Id
            JSONObject obj = view.toJSON();
            obj.put("id", "uuidis_youtube-view");
            view = View.deserialize(obj.toString());
            view.save();
        } catch (Exception ex) {
            ex.printStackTrace();
        }*/
    }
}

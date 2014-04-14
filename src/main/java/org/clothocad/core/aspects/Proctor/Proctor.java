/*
 * 
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
ENHANCEMENTS, OR MODIFICATIONS..
 */

package org.clothocad.core.aspects.Proctor;

import java.util.HashMap;
import java.util.Map;
import java.util.UUID;
import org.clothocad.core.aspects.Aspect;
import org.clothocad.core.aspects.Hopper;
import org.clothocad.core.datums.Doo;
import org.clothocad.core.communication.mind.PageMode;
import org.clothocad.core.datums.ObjectId;
import org.clothocad.model.Person;
import org.clothocad.model.ServerTrailDeprecated;

/**
 * @author John Christopher Anderson
 */


public class Proctor implements Aspect {

    /**
     * SS api-relayed method to start the process of a student taking a
     * quiz to earn a badge
     * 
     * @param student
     * @param quiz
     * @param parentDoo 
     */
    public void initiateTrail(Person student, ServerTrailDeprecated trail, Doo parentDoo) throws Exception {
        //Create a QuizDoo to manage the task
        ProctorDoo doo = new ProctorDoo(parentDoo);
        doo.trailId = trail.getId();
        doo.studentId = student.getId();
        
        //The Trail data contains in the TOC, construct the widget code
        
        //Throw up a new Page with the Trail (via Communicator/Mind)
        //GET THE VIEW FROM THE QUIZ
        //CALL COMMUNICATOR TO ISSUE A SHOW(...) COMMAND FOLLOWED BY AN UPDATE
        //WITH THE DATA HELD AS AN Instance in the quiz
        
        //Put the doo in the hopper
        Hopper.get().add(doo);
    }
//
//    /**
//     * Run the rubric on an answer
//     * 
//     * @param submittedAnswer
//     * @param person
//     * @param doo
//     * @return 
//     */
//    public boolean gradeQuiz(Map<String, Object> submittedAnswer, Person person, Quiz quiz, Doo doo) {
//        ProctorDoo qdoo = new ProctorDoo(doo);
//        qdoo.quizId = quiz.getId();
//        qdoo.studentId = person.getId();
//        
//        /** TODO:
//        Function rubric = quiz.getRubric();
//        try {
//            String resultStr = (String) rubric.execute(submittedAnswer);
//            Map<String, Object> result = new Map<String, Object>(resultStr);
//            
//            //NEED TO SEE WHAT THIS LOOKS LIKE.  WE'RE PROBABLY AT THE POINT WHERE WE NEED TO FIX THIS
//            //YEAH, NEED TO WORK OUT THE SERVERSCRIPT EXECUTION ISSUE
//            
//            
//            
//        } catch (Exception ex) {
//            Logger.getLogger(Proctor.class.getName()).log(Level.SEVERE, null, ex);
//        }
//        **/
//        return false;
//    }
//    
    private class ProctorDoo extends Doo {
        public ProctorDoo(Doo parent) {
            super(parent, true);
        }
        
        ObjectId quizId;
        ObjectId trailId;
        ObjectId studentId;
    }
    
    /**
     * The Proctor can instantiate any kind of page.  Proctor commands are only executed
     * from Trails or as Editor requests.  So, it's under separate wiring than normal views
     * which can only deal with workspaces normally.
     * 
     * @param pageId
     * @param mode
     * @return
     * @throws Exception 
     */
    public Map<String, Object> makeNewPage(PageMode mode) throws Exception {
        Map<String, Object> out = new HashMap<>();
        out.put("mode", mode.toString());
        out.put("ephemeral_link_page_id", UUID.randomUUID().toString());
        out.put("command", "addPage");
        return out;
    }
//
//    /**
//     * Relay method for receiving a show command
//     * @param viewRef
//     * @param sharables
//     * @param sloppy_args 
//     */
//    public final void show(View view, String sharables, String position, Page targetPage){
//        try {
//            
//            //IF THE PROCESS WAS TIED TO A COMMAND-INITIATED PROCESS THERE WOULD BE A PARENT DOO THAT NEEDS TO BE FOUND.
//            Doo parentDoo = Hopper.get().extract(null);
//            ShowDoo doo = new ShowDoo(parentDoo);
//
//            //Gather up referenced objects
//            List<Sharable> shareList = resolveToSharables(sharables);
//            
//            //Create the widget and put it on its page
//            Widget widget = new Widget(targetPage, view);
//            targetPage.addWidget(widget);
//            
//            //Just for record keeping
//            doo.viewId = view.getId();
//            doo.widgetId = widget.getId();
//            doo.collectSharables(shareList);
//            
//            //Create all the commands for the client
//            doo.commandMessageArray = new JSONArray();
//            doo.commandMessageArray.put(   this.makeNewPage(PageMode.TRAILS))
//            doo.commandMessageArray.put(   ServerSideAPI.makeCollect(shareList)           );
//            doo.commandMessageArray.put(   makeShowWidget(widget, position, targetPage.getSocketId()) );
//            doo.commandMessageArray.put(   ServerSideAPI.makeUpdate(widget, shareList)    );
//            doo.commandMessageArray.put(   ServerSideAPI.makeCallback(doo.getId(), targetPage.getSocketId())        );
//            
//            //Put the doo into the hopper to await a callback
//            Hopper.get().add(doo);
//            
//            //Send the commands
//            Communicator.get().sendClientMessage(socket_id,
//                    SendChannels.commandList,
//                    doo.commandMessageArray.toString());
//
//            //Save everything whose state was changed  //JCA:  THIS NEEDS A CALLBACK/FAILURE RESPONSE THAT REVERTS THIS (EVENTUALLY)
//            Persistor.get().persistDatum(mind);
//        } catch (Exception e) {
//            Logger.log(Logger.Level.WARN, "", e);
//            e.printStackTrace();
//            //REVERT THE MIND
//            //SEND CLIENT MESSAGE TELLING THAT IT FAILED
//        }
//    }
    
    /**
     * Instruct the client to display the GUI of a view and position it
     * according to positioning parameters.  This duplicates communicator intentionally.
     * 
     * Since Trails are associated with security in world, the user needs to be confident
     * that when her Clotho shows a Trail that she is watching a Trail and not something
     * being injected from another View.  So, the Proctor handles it's own command objects
     * and there is no means of directing a Widget to a Trails page except by the Proctor.
     * 
     * The position object should be:
     *      int pageIndex; //0, 1, 2, 4, ...  the page number
     *      int x, y, w, h; //and maybe z and depth if ever is 3D
     * 
     * @param viewId
     * @param position
     * @return
     * @throws Exception 
     */
//    private Map<String, Object> makeShowWidget(Widget widget, Map<String, Object> position) throws Exception {
//        Map<String, Object> out = new Map<String, Object>();
//            View view = widget.getView();
//
//            //Load the html and js scripts into a Map<String, Object>, insert the uuid for the widget
//            String widgetId = widget.getId();
//            
//            String html = view.getGraphicsScript(); //this is the correct call that needs to be made in makeShowWidget(...) to get the raw html
//            HTMLParser parser = new HTMLParser(html);  //instantiate some html parser, I've made up this syntax
//            
//            //Iterate through the HTML script, gather up the fieldId's, and replace things with UUIDs
//            for(Element element : parser.getElements()) {
//                 String fieldId = element.getAttribute("id");
//                 //Replace the fieldId in there with the widget's uuid appended to it
//                 element.setAttribute("id", widgetId + fieldId);
//            }
//
//            html = parser.toString();
//            
//            
//            //JUST IGNORE DEALING WITH ONSHOW FOR NOW.
//            String onshow = view.getOnShowScript();
//
//            out.put("widget_id", widgetId);
//            out.put("content", html);
//            out.put("on_show", onshow);
//            out.put("parent_widget_id", "trail_toc");
//            out.put("command", "showWidget");
//            
//        return out;
//    }
    
    //
    //Singleton stuff
    //
    private Proctor() {}
    
    public static Proctor get() {
        return singleton;
    }

    private static final Proctor singleton = new Proctor();
}

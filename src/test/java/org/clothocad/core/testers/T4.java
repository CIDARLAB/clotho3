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
ENHANCEMENTS, OR MODIFICATIONS..
 */

package org.clothocad.core.testers;

import org.clothocad.core.aspects.Collector;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.Instance;
import org.clothocad.core.layers.communication.Callback;
import org.clothocad.core.layers.execution.Executor;
import org.clothocad.core.util.Logger;
import org.json.JSONException;
import org.json.JSONObject;

/**
 * this test pulls the cindy instance of Person, the makeabuddy assistant,
 * then calls makeabuddy on cindy as its sole input.   MakeABuddy creates
 * another person object including the argument's name in their name.
 * So, it will output a new ObjBase representing the buddy
 * 
 * In json, that new object will be something like:
 * {"email":"bob@bobbobcom","name":"cindysus buddy","schemaId":"person_schema_isuuid","uuid":"9c28558e-9f97-421f-bbb2-f04f75557fcd"}
 * 
 * @author John Christopher Anderson
 */
public class T4 {
    public static void main(String[] args) {
        //Retrieve an instance of Person, note that I've manually changed the uuid to
        //something human-readable, but in reality are auto-generated, so requires query-by-name
        Instance Cindy = (Instance) Collector.get().getDatum("specific-cindysu-is-uuid");
        Logger.log(Logger.Level.INFO, "Tester4 has Cindy: " + Cindy.getId());
        
        //Fetch the assistant
        Function asst = (Function) Collector.get().getDatum("specific-makeabuddy-is-uuid");
        if (asst == null) {
            Logger.log(Logger.Level.FATAL, "Could not retrieve 'specific-makeabuddy-is-uuid'");
            return;
        }
        Logger.log(Logger.Level.INFO, "Tester4 has MakeABuddy: " + asst.getDescription());
        
        //Construct the bits needed for Executor start with passing in Cindy
        JSONObject inputs = new JSONObject();
        try {
            inputs.put("cindy", Cindy.toJSON());
        } catch (JSONException e) {
            Logger.log(Logger.Level.FATAL, "JSON error", e);
            return;
        }
        
        System.out.println("tester 4:  inputs: " + inputs.toString());
        
        //Create the asynchronous callback to return data
        Callback callback = new Callback() {
            @Override
            public void onSuccess(JSONObject outputData) {
                System.out.println("Tester4 on success! has " + outputData.toString());
                Logger.log(Logger.Level.INFO, "Yay, it returned successful! with " + outputData.length());
                try {
                    Logger.log(Logger.Level.INFO, outputData.toString());
                } catch(Exception e) {
                    Logger.log(Logger.Level.WARN, "Tester4 result is not a JSONObject", e);
                }
            }

            @Override
            public void onFailure(Throwable err) {
                Logger.log(Logger.Level.WARN, "Bummer, it failed, but didn't crash!");
            }
        };
        
        //Run it
        Executor.get().run(null, asst, inputs, callback);
        Logger.log(Logger.Level.INFO, "I'm finished!");
    }
}

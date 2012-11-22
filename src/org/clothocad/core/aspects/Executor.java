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

package org.clothocad.core.aspects;

import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.Doo;
import org.clothocad.core.datums.util.ServerScript;
import org.json.JSONObject;

/*
 * So, I didn't make up the word Executor, it's an interesting Java interface of thread managers:
 * 
 * http://docs.oracle.com/javase/1.5.0/docs/api/java/util/concurrent/ScheduledThreadPoolExecutor.html
 * http://docs.oracle.com/javase/1.5.0/docs/api/java/util/concurrent/Executor.html
 * http://docs.oracle.com/javase/1.5.0/docs/api/java/util/concurrent/ExecutorService.html
 * http://docs.oracle.com/javase/1.5.0/docs/api/java/util/concurrent/Executors.html
 * 
 * Good description of how to use:
 * http://tutorials.jenkov.com/java-util-concurrent/executorservice.html
 */


/**
 * Executor schedules work (via Process) by Assistants on
 * DataFields (such as ObjBases).
 * @author jcanderson
 */
public final class Executor implements Aspect {
    /**
     * This is the big important method exposed in api to run an assistant
     * on a list of inputs and a callback (asynchronously).
     * 
     * The inputData list is presumably ObjBases, but can also be any of the
     * primitive DataFields like StringField.
     * 
     * The callback duplicates the pattern from GWT of returning one method (onSuccess) carrying the
     *  List<DataField> outputData
     * if successful, and an onFailure with the Exception/Thowable if it has to abort execution.
     * 
     * @param asst  The assistant
     * @param inputData  the value data, usually objbases but can be wrapped primitives and xrefs too
     * @param callback  What to do after the thing has finished (an async callback)
     */
    public synchronized final void run(Doo parent, final Function assistant, JSONObject jsonInputs, Callback callback) {
        ExecutorDoo doo = new ExecutorDoo(parent, assistant, jsonInputs, callback);
        doo.setMessage("Executor Starting to run doo");
        if(runCanDooIt(doo)) {
            doo.setMessage("Executor says it can runCanDooIt");
            JSONObject result = runDooit(doo);
            if(result==null) {
                doo.setMessage("Executor results of runDooIt were null");
                callback.onFailure(null);
            } else {
                doo.setMessage("Executor results of runDooIt were successful");
                callback.onSuccess(result);
            } 
            doo.terminate();
        }
    }
    
    //7/1 correct
    private synchronized boolean runCanDooIt( ExecutorDoo doo ) {
        try {
            //Get the language and script
            doo.setMessage("Executor.runCanDooIt about to test");
            ServerScript candooit = doo.assistant.getCanDooIt();
            
            
            String resultStr = candooit.run(doo.inputs);
            Boolean bool = Boolean.parseBoolean(resultStr);
            return bool;
        } catch (Exception e) {
            doo.setMessage("Executor.runCanDooIt threw an exception for some reason");
            Logger.log(Logger.Level.WARN,
                       "This needs to relay a developer message as the candooit method failed to run properly",
                       e);
            doo.terminate();
            return false;
        }
    }
    
    private synchronized String escape(String input) {
        StringBuilder sb = new StringBuilder();
        for(int i=0; i<input.length(); i++) {
            char achar = input.charAt(i);
            if(achar=='"') {
                sb.append('\\');
            }
            sb.append(achar);
        }
        return sb.toString();
    }
    
    //7/1 correct
    private synchronized JSONObject runDooit(ExecutorDoo doo) {
        try {
            //Get the language and script
            ServerScript dooit = doo.assistant.getDooIt();
            Logger.log(Logger.Level.INFO, "rundooit about to call run");
            String resultStr = dooit.run(doo.inputs);
            JSONObject result = new JSONObject(resultStr);
            return result;
        } catch (Exception e) {
            doo.setMessage("Executor.runCanDooIt threw an exception for some reason");
            Logger.log(Logger.Level.WARN,
                       "This needs to relay a developer message as the rundooit method failed to run properly",
                       e);
            doo.terminate();
            return null;
        }
    }

    /**
     * An instance of Doo is instantiated on a per command execution task.  It holds
     * state for the elements of the execution and manages its lifecycle.  Eventually
     * Executor will have another method for scheduling Doo's to a specific datetime
     */
    public class ExecutorDoo extends Doo {
        final Function assistant;
        final JSONObject inputs;
        final Callback callback;

        public ExecutorDoo(Doo doo, Function assistant, JSONObject inputs, Callback callback) {
            super(doo, false);
            this.assistant = assistant;
            this.inputs = inputs;
            this.callback = callback;
        }
    }
    
    public static Executor get() {
        return singleton;
    }
    private static final Executor singleton = new Executor();
}

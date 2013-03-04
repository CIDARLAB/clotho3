package org.clothocad.app.GenomeRefactoring;

import java.util.ArrayList;
import java.util.List;

import org.clothocad.core.aspects.Executor;
import org.clothocad.core.aspects.Persistor;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.Schema;
import org.clothocad.core.datums.objbases.Person;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.datums.util.FieldType;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.datums.util.ServerScript;
import org.clothocad.core.layers.communication.Callback;
import org.clothocad.core.util.FileUtils;
import org.clothocad.core.util.Logger;
import org.json.JSONException;
import org.json.JSONObject;

public class GenomeRefactoring {

	public GenomeRefactoring() {
		testFunctions();		
		// testScripting();
	}
	
	private void testFunctions() {
		
		executeFunction(
				createScript());
		
		// first, create the functions
		//ArrayList<Function> lstFunctions = createFunctions();

		// now, execute the functions
		//Function f1 = lstFunctions.get(0);
		//this.executeFunction(f1);
	}
	
	private Function createScript() {
		
        ServerScript dooit = new ServerScript(
        		FileUtils.readFile("./scripts/genome-refactoring.py"),
                Language.python);
		/**
        ServerScript dooit = new ServerScript(
        		FileUtils.readFile("./scripts/sequence_transformation.py"),
                Language.python);
      	**/
		
        // can doo it
        ServerScript candooit = new ServerScript("return true;", Language.JavaScript);
        
        // schema
        ///Schema personSchema = Person.getSchema();
        
        // input args
        List<ClothoField> inputArgs = new ArrayList<ClothoField>();
        ClothoField uuidSeq = new ClothoField(
        		"sequence", 
        		FieldType.STRING, 
        		"ATCG", 
        		1);
        inputArgs.add(uuidSeq);
            
        // output args
        List<ClothoField> outputArgs = new ArrayList<ClothoField>();
        /***
        ClothoField reverse = new ClothoField(
        		"sequence", 
        		FieldType.STRING, 
        		personSchema.getId(), 
        		1);
        outputArgs.add(reverse);
		***/
        
        return Function.create(
        		Person.getAdmin(), 
        		"remove_feature_overlap", 
        		"this is the remove_feature_overlap function's documentation", 
        		dooit, 
        		candooit, 
        		inputArgs, 
        		outputArgs);
	}
	
	private boolean executeFunction(Function f) {
		// 1. define the function's inputs
        JSONObject sequence = new JSONObject();
        try {
            sequence.put("sequence","ATCAAAAAGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAATTACCCCCCAAGTCTA");
        } catch (JSONException e) {
            Logger.log(Logger.Level.FATAL, "JSON error", e);
            return false;
        }
        
        // 2. define a callback that process the function's response
        Callback callback = new Callback() {
            @Override
            public void onSuccess(JSONObject outputData) {
                Logger.log(Logger.Level.INFO, "Yay, it returned successful! with " + outputData.length());
                try {
                    Logger.log(Logger.Level.INFO, outputData.toString());
                } catch(Exception e) {
                    Logger.log(Logger.Level.WARN, "remove_feature_overlaps result is not a JSONObject", e);
                }
            }

            @Override
            public void onFailure(Throwable err) {
                Logger.log(Logger.Level.WARN, "Bummer, it failed, but didn't crash!");
            }
        };
        
        // 3. execute the function
        Executor.get().run(null, f, sequence, callback);
        Logger.log(Logger.Level.INFO, "I'm finished!");
        return true;
	}

	private void testScripting() {
		// just execute a python script		
		ServerScript script = new ServerScript(
				FileUtils.readFile("./scripts/genome-refactoring.py"), 
				Language.python);
		
		try {
			script.run(null);
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
	
	public static void main(String[] args) {
		new GenomeRefactoring();
	}

}

package org.clothocad.app.GenomeRefactoring;

import java.util.ArrayList;
import java.util.List;

import org.clothocad.core.aspects.Collector;
import org.clothocad.core.aspects.Executor;
import org.clothocad.core.aspects.Persistor;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.Instance;
import org.clothocad.core.datums.Schema;
import org.clothocad.core.datums.objbases.Person;
import org.clothocad.core.datums.util.ClothoField;
import org.clothocad.core.datums.util.FieldType;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.datums.util.ServerScript;
import org.clothocad.core.layers.communication.Callback;
import org.clothocad.core.util.Logger;
import org.json.JSONException;
import org.json.JSONObject;

public class ReverseComplement {

	public ReverseComplement() {
		//createReverseComplementFunction();
		testReverseComplementFunction();
	}
	
	private void createReverseComplementFunction() {		
        // I'm an App developer and want to create a new Function 
		// that produces the reverse complement of a DNA sequence
		// -> therefore I have to describe the functions behavior

		// the function takes as input argument a UUID
		// and returns the UUIDs reverse complement (if the UUID exists and the UUID contains a sequence)
		
		// we have to create a ServerScript Doo
        ServerScript dooit = new ServerScript(
                "clotho.log(\"INFO\", '+++running Dooit in js from reverseComplement');\n"
              + "clotho.log(\"INFO\", '+++calculating reverse complement of: ' + uuid);\n"
              + "return uuid;",
                 Language.JavaScript);
      
        // can doo it
        ServerScript candooit = new ServerScript("return true;", Language.JavaScript);
        
        // schema
        Schema personSchema = Person.getSchema();
        
        // input args
        List<ClothoField> inputArgs = new ArrayList<ClothoField>();
        ClothoField cindy = new ClothoField("uuid", FieldType.STRING, personSchema.getId(), 1);
        inputArgs.add(cindy);
            
        // output args
        List<ClothoField> outputArgs = new ArrayList<ClothoField>();
        ClothoField reverse = new ClothoField("reverse", FieldType.STRING, personSchema.getId(), 1);
        outputArgs.add(reverse);

        Function fReverseComplement = Function.create(
        		Person.getAdmin(), 
        		"reverseComplement", 
        		"this function returns the reverse complement of a given UUID", 
        		dooit, 
        		candooit, 
        		inputArgs, 
        		outputArgs);
        
        // now, store the reverseComplement() function in the DB
        Persistor.get().persistDatum(fReverseComplement);
	}
	
	private void testReverseComplementFunction() {
		// now, I'm a Clotho User and 
		// want to utilize to reverseComplement() function

		// first, we load a specific DNA component from the DB, i.e. GFP
        Instance gfp = (Instance) Collector.get().getDatum("specific-gfpuv-is-uuid");
        try {
			Logger.log(Logger.Level.INFO, "SimpleFeature: " + gfp.getString("sequence"));
		} catch (JSONException e) {
			e.printStackTrace();
			System.exit(1);
		}
		
        /* Question:
         * How can we retrieve the reverseComplement() function if the 
         * function's UUID was created automatically? 
         */
        
        // then, we fetch the assistant
        Function fReverseComplement = (Function) Collector.get().getDatum("64e27a0d-888c-4512-a3a2-bfde8e760092");
        if (fReverseComplement == null) {
            Logger.log(Logger.Level.FATAL, "Could not retrieve '64e27a0d-888c-4512-a3a2-bfde8e760092'");
            return;
        }
        Logger.log(Logger.Level.INFO, "reverse complement function description: " + fReverseComplement.getDescription());
        
        // now, we have to create the function's input arguments 
        JSONObject inputs = new JSONObject();
        try {
            inputs.put("uuid", gfp.toJSON());
        } catch (JSONException e) {
            Logger.log(Logger.Level.FATAL, "JSON error", e);
            return;
        }
        
        System.out.println(fReverseComplement.getId()+" inputs: " + inputs.toString());
        //Executor.get().run(null, fReverseComplement, inputs, null);

        /***
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
        ***/
        
        Logger.log(Logger.Level.INFO, "I'm finished!");
	}
	
	public static void main(String[] args) {
		new ReverseComplement();
	}

}

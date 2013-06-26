package org.clothocad.core.testers.api;

import java.util.Random;

import org.clothocad.core.layers.communication.Channel;
import org.clothocad.core.layers.communication.Router;
//import org.clothocad.dm.clotho2.NucSeq;
import org.json.JSONObject;

public class Creator {

    /***
	private static final int NR_OF_SEQUENCES = 1000;
	
	public Creator() {
		
	}
	
	public static void main(String[] args) 
			throws Exception {
	
		Creator creator = new Creator();
		
		// first, we create some JSON objects 
		// and insert them into the database
		for(int i=1; i<=NR_OF_SEQUENCES; i++) {
			creator.insertNucSeq();
		}
	}
	
	public void insertNucSeq() 
			throws Exception {		
		NucSeq ns = new NucSeq(this.generateRandomSequence());
		
		// NucSeq.toJSON();
		
		// now, we need to create a JSON object 
		JSONObject data = new JSONObject();
		try {
			data.put("uuid", ns.getUUID().toString());
			data.put("className", NucSeq.class.getCanonicalName());
			JSONObject model = new JSONObject();
			model.put("sequence", ns.getSeq());
			data.put("model", model);

			// now we say Clotho.set();
			// but first, we send the NucSeq JSON object to the Router directly
			Router.get().receiveMessage(null, Channel.create.toString(), data);				
		} catch(Exception e) {
			throw new Exception(e);
		}
	}
	
	public String generateRandomSequence() {
		
		StringBuilder sb = new StringBuilder();
		Random rand = new Random();
		for(int i=1; i<=1000; i++) {
			int n = rand.nextInt()%4;
			switch(n) {
			case 0:
				sb.append("A");
			case 1:
				sb.append("T");
			case 2:
				sb.append("C");
			case 3:
				sb.append("G");
			}
		}
		return sb.toString();
	}
  ***/
}

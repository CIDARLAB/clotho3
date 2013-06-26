package org.clothocad.core.testers.api;

import java.util.ArrayList;
import java.util.List;

import org.clothocad.core.layers.communication.Channel;
import org.clothocad.core.layers.communication.Router;
//import org.clothocad.dm.clotho2.NucSeq;
import org.json.JSONException;
import org.json.JSONObject;

public class Getter {
/***
	public static void main(String[] args) {
		Getter g = new Getter();
		
		// get all nucseq objects
		List<NucSeq> lst = g.getAllNucSeqs();
		
		g.getNucSeqById(lst);
	}
	
	public List<NucSeq> getAllNucSeqs() {

		JSONObject data = new JSONObject();

		try {
			data.put("className", NucSeq.class.getCanonicalName());

			// now we say Clotho.set();
			// but first, we send the NucSeq JSON object to the Router directly
			Router.get().receiveMessage(null, Channel.get.toString(), data);				
		} catch (JSONException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return new ArrayList<NucSeq>();

	}
	
	public void getNucSeqById(List<NucSeq> lst) {
		if(null != lst) {
			// iterate over all nucseqs in the list
			for(NucSeq ns : lst) {
				System.out.println(ns);
			}
		}
	}
	
	public NucSeq getById(String uuid) {
		JSONObject data = new JSONObject();
		try {
			// either:
			data.put("uuid", uuid);
			// now we say Clotho.set();
			// but first, we send the NucSeq JSON object to the Router directly
			Router.get().receiveMessage(null, Channel.get.toString(), data);
			
			return new NucSeq();
			
		} catch(Exception e) {
			e.printStackTrace();
		}
		return null;
	}
***/ 
}

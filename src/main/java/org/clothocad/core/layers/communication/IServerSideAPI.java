package org.clothocad.core.layers.communication;

import org.clothocad.core.datums.Sharable;
import org.clothocad.core.datums.View;
import org.json.JSONObject;

public interface IServerSideAPI {

	public JSONObject get(String uuid);

	public void set(String sSharableID, JSONObject objValues);
	public void create(Class cType, JSONObject objValues);
	
	public void edit(); // -> what's the difference to set() ??
	
	public void revert(); // unclear
	
	public void say(String sMessage); // who says what to whom? 
	
	public void alert(String sMessage); // 
	
	public void show(Sharable objSharable, View objView);
	
	public void run();
	
	public void newPage();
	public void closePage();
	
	public void login();
	public void logout();
	public void lern();  // for the learning trails -> JCA	

	// missing functionalities:
	// - create() .. to insert new datums/sharables into the db
	// or should this be done in the set method?
	// e.g. if the sSharableID does not exist yet, insert it into the db...
	
	// - delete/remove()
}

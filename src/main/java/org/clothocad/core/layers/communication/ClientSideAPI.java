package org.clothocad.core.layers.communication;

import org.clothocad.core.datums.ObjBase;

public interface ClientSideAPI {

	
	 /*
	  * This is part of the middleware-communications between clotho and its clients.  
	  * It is not accessible to the ssAPI.  The ssAPI is invoked universally--REST, search bar, server scripts, 
	  * view scripts, etc. all speak to Clotho through the same API.  So, dynamic content is never 
	  * directly aware of these commands.  
	  * 
	  * When a client does a get request, the object would be sent back to the client's collector via a collect call.  
	  * So, it isn't invoked directly, but it happens.  
	  * This is for caching of the object.
	  *   
	  * So, the process should ideally go:  
	  * 1. server sends the JSON to the client as a collect call, 
	  * 2. client puts the JSON into the collector, 
	  * 3. client returns that JSON to the calling function (sync) or to the callback of the calling function (async).
	  */
	 public void collect(ObjBase obj);
	 
	 
	 /*
	  * This *may* be a client-side api method.  It's not an ssAPI method.  I'm not sure it is required at all, though.  
	  * It was intended to trigger the broadcast of the data change to widgets when we weren't using a real MVC.  
	  * Now that GUI changes will be automatically bound to data, only the 'collect' method should be required for 
	  * MVC auto-updates. 
	  * 
	  * The data just changes due to a push of JSON. 
	  * 
	  * The data listeners respond automatically.
	  */
	 public void update(ObjBase obj);

	 /*
	  * showWidget
	  *
	  * 	 see above discussion of show(...).

	  */


	 /*
	  * remove
	  *
	  * I'm not sure if this means "remove this widget from the page" which used to be called removeWidget in librecv.  
	  * That's a clientside API method.  The other 'remove' might be 'delete' meaning 
	  * "check all the dependencies for this sharable, delete the ones you won't need anymore, then delete that sharable". 
	  * This delete call would be an ssAPI method, and it can be refused by permissions like many of these ssAPI calls.
	  */
}


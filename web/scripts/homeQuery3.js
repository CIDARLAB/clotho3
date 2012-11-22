
/*
//given a valid parsed JSON, adds results and relevant suggestions to the main dialog
function handle_submit(json) {
	$('#main-container').removeClass("ajax-loading");
	message('<div class="dialog_update">Sent: '+json.query);
	
	message('<div class="dialog_response">Command: ' + json.result.command);
	message('<div class="dialog_response">Result: ' + json.result.result);
	
	if (json.suggestions.length > 0) {
		$('#main-content-dialog').append('<div class="dialog_suggestionList"><div class="suggestionToggler"></div><div class="suggestions_container"></div></div>');
		$.each(json.suggestions, function(i, object) {
			$('.suggestions_container:last').append('<a class="dialog_suggestion">'+json.suggestions[i].command+'</a>');
		});
	}
	
	$('.dialog_suggestionList:last .suggestions_container').hide();
	$('.dialog_suggestionList:last a.dialog_suggestion').click(function(event) {
		event.preventDefault();
		//fix to grab content between <a></a>
		var command = $(this).text();
        var message = JSON.stringify({"command":command, "query":getQuery});
        libsend.call("submitCommand", message);
	});
	$(".dialog_suggestionList:last .suggestionToggler").click(function(){
		$(this).next(".suggestions_container:last").slideToggle(400);
	});
}
 
// **REQ** html_string - html to be shown
// **REQ** uuid - not view uuid, uuid generated and passed from server 
//pos - position to show on screen - set as element.style
//mode - how the object should be shown:
//	0 = Homepage
//	1 = Trails
//	2 = BrowseMode
//	3 = Editor
//	4 = Workspace
function show(html_string, uuid, posx, poxy, objWidth, objHeight, mode) {
	//defaults for variables
	//NOTES
	//not dealing with more right now
	//do we want better positioning logic? i.e. take size of widget into account
	posx = typeof posx !== 'undefined' ? posx : 20;
	posy = typeof posy !== 'undefined' ? posy : 20;
	objWidth = typeof objWidth !== 'undefined' ? objWidth : 200;
	objHeight = typeof objHeight !== 'undefined' ? objHeight : 230;
	
	//THIS IS A PLACEHOLDER IF UUID AND HTML_STRING ARE NOT PASSED
	//LATER THIS SHOULD BE ERROR HANDLING
	uuid = typeof uuid !== 'undefined' ? uuid : Math.floor(10000000*Math.random());
	html_string = typeof html_string !== 'undefined' ? html_string : '<div class="default-form"> \n\
				<div class="text-div name-div"> \n\
					<label for=""><strong class="name-label">Name</strong></label> \n\
					<input type="text" name="Name" class="Name" /> \n\
				</div> \n\
				<div class="text-div email-div"> \n\
					<label for="Email"><strong class="email-label">Email</strong></label> \n\
					<input type="text" spellcheck="false"  name="Email" class="Email" value="" /> \n\
				</div> \n\
				<input type="submit" class="blue_button update" name="update" value="Update" onclick=\'update('+uuid+')\' /> \n\
			</div>';
	
	//assign as element.style so that can be changed by jquery when interacted with (e.g. drag, resize)
	//default z-index to 1000 so shows up on top. will be reset to appropriate lower number when dragged
	var widget = $(
		'<div uuid="' + uuid + '" class="widget resizable" style="z-index: 1000; position: absolute; top: '+posy+'; left: '+posx+'; width: '+objWidth+'px; height: '+objHeight+'px;"> \n\
			<div class="widget_close_button"></div> \n\
			<div class="widget-container"> \n\
				' + html_string + ' \n\
			</div> \n\
		</div>');
		
	$('#widget_workspace').append(widget);
	activateWidget(widget);
}

function hide(uuid) {
	//removes all jQuery handlers associated with the element
	$("div[uuid="+uuid+"]").remove();
}

//grab widgetsin jquery with CSS selector like this: $("div[uuid=<unique_uuid>]")
//This line sets the content of the widget to newData... here just using the example JSON to populate it
//$("div[uuid="+uuid+"] .widget-container").html(newData);
function update(uuid, newData) {
	//need to make this more dynamic
	// how do we know what data is being passed? Is this just a basic method that is used to interact with the server, but each widget will have its own specialized update method?
	
	//THIS IS A PLACEHOLDER IF UUID AND HTML_STRING ARE NOT PASSED
	//LATER THIS SHOULD BE ERROR HANDLING
	newData = typeof newData !== 'undefined' ? newData : '{"result": {"name": "'+uuid+'", "email": "this@tester.com"} }';
	
	var update_json = JSON.parse(newData);
		
	$("div[uuid="+uuid+"] .Name").val(update_json.result.name);
	$("div[uuid="+uuid+"] .Email").val(update_json.result.email);
}
*/

//need to figure out what to do here... previously was dumping old content and replacing with an AJAX loading icon
//may want to store old content in a variable, dump it, and bring it back if query is empty and not focused on command bar
function presend() {}


/***************
Command Bar Helpers
***************/
	
/*
//appends to dialog
function message(msg){  
	$('#main-content-dialog').append(msg+'</div>');  
}
function message_a(msg){  
	$('#main-content-dialog').append(msg+'</a>');  
}
function message_p(msg){  
	$('#main-content-dialog').append(msg+'</p>');  
}

//append to suggestions
function suggest(msg){  
	$('#theQuery-suggestions ul').append(msg+'</li>');  
}
function suggest_p(msg){  
	$('#theQuery-suggestions').append(msg+'</p>');  
}
*/
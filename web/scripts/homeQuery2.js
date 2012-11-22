// AJAX version
$(document).ready(function(){

	//default session_id
	var session_id = 123456789;
	
	//sends full string, returns results + suggestions to main dialog. called by button and pressing enter
	function send(suppliedQuery){    
		try{  
			if (!suppliedQuery) {
				var curQuery = getQuery();
			} else {
				curQuery = suppliedQuery;
			}
			
			$.ajax({
				//change this to something relevant
				url: "http://www.audiobooks.org/JSON-return.php",
				cache: false,
				type: "get",
				data: jsonify_send(suggestions, curQuery, session_id),
				beforeSend: function() {
					$('#main-content-default').fadeOut("fast");
					$("#main-container").addClass("ajax-loading");
				},
				success: function(response, textStatus, jqXHR){
					$('#main-container').removeClass("ajax-loading");
					message('<div class="dialog_update">Sent: '+curQuery);
					
					//handle the passed JSON here... needs to be valid JSON to display
					var passed_json = jQuery.parseJSON(response);
					var num_suggestions = passed_json.message.suggestions.length
					//
					//make this not each() later (unless pass multiple results...)
					//
					$.each(passed_json.messsage, function(i, object) {
						message('<div class="dialog_response last">Command: '+object.result);
						if (num_suggestions > 0) {
							$('#main-content-dialog').append('<div class="dialog_suggestionList">');
							$.each(object.suggestions, function (i, object) {
								message('<a class="dialog_suggestion">'+object.command);
							});
							$('#main-content-dialog').append();
						}
					});
				},
				//make suggestions executable
				complete: function(jqXHR, textStatus) {
					$('a.dialog_suggestion').click(function(event) {
						event.preventDefault();
                        //fix to grab content between <a></a>
						var identity = $(this).attr('id');
						send(identity);
					});
				},
				error: function(jqXHR, textStatus, errorThrown){
					message_p('<p class="error">AJAX Error: ' + errorThrown + textStatus); 
				}
			});

		} catch(exception){  
			message_p('<p class="error">General Error: ' + exception);  
		}  
		//$('#theQuery').val("");
	} 
	
	//get recommended queries for passed query fragments, called on keyup, results pushed to suggestions bar below search bar
	function get_suggestions(suppliedQuery){
		try{  
			if (!suppliedQuery) {
				var curQuery = getQuery();
			} else {
				curQuery = suppliedQuery;
			}
			
			$.ajax({
				//change this to something relevant
				url: "http://www.audiobooks.org/JSON-return.php",
				cache: false,
				type: "get",
				data: jsonify_send(suggestions, curQuery, session_id),
				beforeSend: function() {
				},
				success: function(response, textStatus, jqXHR){
					var passed_json = jQuery.parseJSON(response);
					var num_results = passed_json.result.length;
					$("#theQuery-suggestions").html('<ul></ul>').show();
					$.each(passed_json.message, function(i, object) {
						suggest('<li'+(i == num_results-1 ? ' class="suggestion last"' : '')+'>'+object.command+'<div class="suggestion_confidence">'+(object.confidence ? object.confidence+'%' : '')+'<div>');
					})
				},
				//make suggestions executable
				complete: function(jqXHR, textStatus) {
					$('a.theQuery-command').click(function(event) {
						event.preventDefault();
						var identity = $(this).attr('id');
						send(identity);
					});
				},
				error: function(jqXHR, textStatus, errorThrown){
					$("#theQuery-suggestions").show();
					suggest_p('<p class="error">AJAX Error: ' + errorThrown + textStatus); 
				}
			});

		} catch(exception){  
			$("#theQuery-suggestions").show();
			suggest_p('<p class="error">General Error: ' + exception);  
		} 
	}
	
	function jsonify_send(channel, message, session){
		return "{channel: "+channel+", message: "+message+", session: "+session+"}";
	}
	
	//appends div to dialog
	function message(msg){  
		$('#main-content-dialog').append(msg+'</div>');  
	}
	//appends paragraph to dialog
	function message_p(msg){  
		$('#main-content-dialog').append(msg+'</p>');  
	}
	
	//append suggestion as list element
	function suggest(msg){  
		$('#theQuery-suggestions ul').append(msg+'</li>');  
	}
	//replace suggestions with paragraph
	function suggest_p(msg){  
		$('#theQuery-suggestions').append(msg+'</p>');  
	}
	
	function getQuery() {
		//later need to escape this
		var query = $("#theQuery").val();
		if (query == "Enter your command or search term") {
		   query = "DEFAULT PHRASE";
		}
		if (query == "") {
		   query = "EMPTY";
		}
		return query;
	}
	
	//
	// Key-bindings
	//
	
	//prevent leaving page when press enter
	$('#theQuery').keypress(function(event) {  
		if (event.keyCode == '13') {  
			event.preventDefault();
			send();
		}  
	});
	
	//get suggestions after each keystroke
	$("#theQuery").keyup(function () {
		get_suggestions();
	});
	
	
	//clicking search button launches seach also
	$("#theQuery-button").click(function() {
		send();
	});
	
	//mousing out of command-bar hides suggestions
	$("#theQuery").blur(function() {
		$('#theQuery-suggestions').fadeOut("fast");
	});
	
	$("#theQuery").focus(function() {
		$('#theQuery-suggestions').fadeIn("fast");
	});
});
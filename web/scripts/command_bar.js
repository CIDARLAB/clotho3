$(window).unload(function () {
    libsend.call("unlinkPage", "");
});

$(document).ready(function(){
	//want to call this to activate the examples
	//remove this once the server passes them
	command_bar.activateSuggestions();
	
	/***************
	Key Bindings
	***************/
	
	//prevent leaving page when press enter
	$('#theQuery').keypress(function(event) {  
		if (event.keyCode == '13') {  
			event.preventDefault();
			$('#theQuery-suggestions').fadeOut("fast");
			command_bar.submitCommand();
		}  
	});
	
	//get suggestions after each keystroke
	$("#theQuery").keyup(function () {
        libsend.call("updateQuery", command_bar.getQuery());
	});
	
	
	//clicking search button launches seach also
	$("#theQuery-button").click(function (e) {
		e.preventDefault();
		command_bar.submitCommand();
	});
	
	//mousing out of command-bar hides suggestions
	$("#theQuery").blur(function() {
		$('#theQuery-suggestions').fadeOut("fast");
		clearTimeout(suggHideTimer);
	});
	
	//show suggestions when focus on command bar
	$("#theQuery").focus(function() {
		$('#theQuery-suggestions').fadeIn("fast");
		clearTimeout(suggHideTimer);
	});
	
	//set binding so mousing outside of suggestions will close it after a delay
	var suggHideTimer;
	$("#command-bar-container").on("mouseleave", (function(e) {
		suggHideTimer = setTimeout('$("#theQuery-suggestions").slideUp(400); $("#theQuery").blur()', 800);
	}));
	
	$("#theQuery-suggestions").on("mouseenter", (function(e) {
		clearTimeout(suggHideTimer);
	}));
});

var command_bar = new Object();

//gets current query
command_bar.getQuery = function() {
	//later need to escape this
	var query = $("#theQuery").val();
	if (query == "Enter your command or search term") {
	   query = "DEFAULT PHRASE";
	}
	return query;
};

command_bar.setQuery = function(newQuery) {
	$("#theQuery").val(newQuery);
};

command_bar.submitCommand = function(command) {
	var new_command = typeof command === 'undefined' ? command_bar.getQuery() : command;
	var message = JSON.stringify({"command":new_command, "query":command_bar.getQuery()});
	libsend.call("submitCommand", message);
	command_bar.postArticle(command_bar.wrap_dialog_bubble(command_bar.wrap_dialog_block('Posted to server: ' + new_command), "dialog_post", "rightArrow"));
	command_bar.scrollDialog();

};

command_bar.showSuggestions = function(suggestions) {
    $("#theQuery-suggestions").html('<ul></ul>').show();
    for (i in suggestions) {
        //could use CSS selectors to simplify this... i.e. li.suggestion:last-child
        $('#theQuery-suggestions ul').append('<li class="suggestion"><a class="theQuery-command">' + suggestions[i].command + '</a></li>');
    }
   
   command_bar.activateSuggestions();
};

//separate for example if we want to activate the home page ones for example
command_bar.activateSuggestions = function() {
	 $('a.theQuery-command').click(function(event) {
        event.preventDefault();
		//submit
		if ($(this).text() != $('#theQuery').val()) {
			command_bar.submitCommand($(this).text());
		}
		
		// activate to set query text to selected command 
			//suggestions at present don't change, just hide.... but set the initial text to the query. if user makes a keystroke, however, suggestions will change
			//this may throw off learning. only the first query will show the partially complete text entered by the user if the line below is uncommented
		command_bar.setQuery($(this).text());
    });
	
	
}

command_bar.scrollDialog = function() {
	// max-height set in CSS
	if ($('#main-content-dialog').height() >= 400) {
		
		//this way needs to be fixed if we want to calculate this
		//var wantedTop = $('.dialog_response:last-child').position().top;
		
		var dialogHeight = $('#main-content-dialog')[0].scrollHeight;
		//scroll to bottom
		$('#main-content-dialog').animate({ scrollTop: dialogHeight }, "slow", 'easeOutBounce')
	}
}

//posts an article... requires that bubbles be inside
//use this because can't post partial divs, so this  makes wrapping easier
command_bar.postArticle = function(message) {
	var post = '<div class="dialog_article">' + message + '</div>';
	$('#main-content-dialog').append(post);
}

//wraps text in a bubble.... does not post
//second two paramters are optional
//side is the side with the arrow - Use the class e.g. rightArrow or leftArrow etc.
//additional classes should be separated by spaces (e.g. "dialog_post other_class")
command_bar.wrap_dialog_bubble = function(message, additionalClasses, side) {
	var arrowText = typeof side !== 'undefined' ? '<span class="arrow ' + side + '"></span>' : '';
	var additionalClasses = typeof additionalClasses !== 'undefined' ? ' ' + additionalClasses : '';
	return '<div class="dialog_bubble' + additionalClasses + '">' + arrowText + message + '</div>';
}

//wraps message in dialog_response - generally will want to use wrap_dialog_bubble above and supply extra class
command_bar.wrap_dialog_response = function(message) {return '<div class="dialog_response">' + message + '</div>';}

//wraps message in dialog_post - generally will want to use wrap_dialog_bubble above and supply extra class
command_bar.wrap_dialog_post = function(message) {return '<div class="dialog_post">' + message + '</div>';}

//wrap in dialog_block
command_bar.wrap_dialog_block = function(message) {return '<div class="dialog_block">' + message + '</div>';}
 $(document).ready(function(){
   //check that browser supports websockets
	if(!("WebSocket" in window)){  
		$("#theQuery-button").click(function() {
			$('#main-content').fadeOut("fast");  
			$('<p>Oh no, you need a browser that supports WebSockets. How about <a href="http://www.google.com/chrome">Google Chrome</a>?</p>').appendTo('#main-container'); 
		});
	} else
	{
		connect();

		function connect(){
			var socket;
			var host = "ws://localhost:8000/socket/server/daemon_start.php";
			
			try {
				socket = new WebSocket(host);
				message('<p class="event">Socket Status: '+socket.readyState);  
	
				socket.onopen = function(){  
					 message('<p class="event">Socket Status: '+socket.readyState+' (open)');  
				}; 
				socket.onmessage = function(msg){  
					 message('<p class="message">Received: '+msg.data);
					 //built-in method handles JSON objects.... e.g. {"name":"John"} can be called as: passedobj.name == "John" will return true
					 var passedobj = jQuery.parseJSON(msg);
				}; 
				socket.onclose = function(){  
					 message('<p class="event">Socket Status: '+socket.readyState+' (Closed)');  
				};
			} catch(exception) {
				message('<p>Error: '+exception);
			}
			
			//sends query to the websocket.
			//default is no argument - gets what is in commandBar
			//can supply query as argument e.g. for running example
			function send(suppliedQuery){    
				try{  
					if (!suppliedQuery) {
						var curQuery = getQuery();
					} else {
						curQuery = suppliedQuery;
					}
					//don't send if empty
					if (curQuery != "EMPTY") {
						$('#main-content-default').fadeOut("fast");
						$("#main-container").addClass("ajax-loading");
						socket.send(curQuery);
						$('#main-container').removeClass("ajax-loading");
						message('<p class="event">Sent: '+curQuery);
					}
				} catch(exception){  
					message('<p class="warning"> Error:' + exception);  
				}  
				//$('#theQuery').val("");
			}
			
			//posts message to appropriate area and adds paragraph close tag (so omit when call this function)
			function message(msg){  
				$('#main-content-dialog').append(msg+'</p>');  
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
			
			/***********
			key and event bindings
			***********/
			
			//prevent leaving page when press enter
			$('#theQuery').keypress(function(event) {  
				if (event.keyCode == '13') {  
					event.preventDefault();
				}  
			});
			
			//send after each keystroke
			$("#theQuery").keyup(function () {
				send();
   			});
			
			//clicking search button launches seach also
			$("#theQuery-button").click(function() {
				send();
			});
			
		} 
	}
});
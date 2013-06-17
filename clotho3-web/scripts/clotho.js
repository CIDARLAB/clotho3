			var server = {
				connect : function() {

					//var location = "ws://localhost:61623";
					var location = "ws://localhost:8080/websocket";

					this._ws = new WebSocket(location);
					this._ws.onopen = this._onopen;
					this._ws.onmessage = this._onmessage;
					this._ws.onclose = this._onclose;
				},
		
				_onopen : function() {
					console.log('SOCKET OPENED');
				},
		
				_send : function(channel, data) {
					
					if (this._ws) {				        
						var message = { 
							'channel': channel, 
							'correlation': '1234567890', 
							'authentication': 'XXXXXX', 
							'data': data 
						};
						
						console.log(message);
						
						this._ws.send(
								JSON.stringify(message));
					}
				},
		
				send : function(channel, message) {
					server._send(channel, message);
				}, 

				_onmessage : function(m) {
					if (m.data) {
						alert(m.data);
					}
				},
		
				_onclose : function(m) {
					this._ws = null;
				}
			};

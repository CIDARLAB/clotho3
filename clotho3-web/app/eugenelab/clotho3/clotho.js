			var server = {
				connect : function() {

					var location = "ws://localhost:8080/websocket";

					this._ws = new WebSocket(location);
					this._ws.onopen = this._onopen;
					this._ws.onmessage = this._onmessage;
					this._ws.onclose = this._onclose;
				},
		
				_onopen : function() {
					console.log('SOCKET OPENED');
				},
		
				_send : function(channel, action, data) {
					
					if (this._ws) {				        
						var message = { 
							'channel': channel, 
							'action': action, 
							'correlation': '1234567890', 
							'authentication': 'XXXXXX', 
							'data': data 
						};
						
						console.log(message);
						
						this._ws.send(
								JSON.stringify(message));
					}
				},
		
				send : function(channel, action, message) {
					server._send(channel, action, message);
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

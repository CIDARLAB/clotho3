/**
 * @description
 * $clotho.socket is the actual WebSocket
 * $clotho.$socket is the socket service
 */

angular.module('clotho.core').service('Socket',
	function($window, $q, $log, PubSub, ClientAPI, Debug) {

	//note - ensuring page-wide singleton
	return ($window.$clotho.$socket) ?
			$window.$clotho.$socket :
			$window.$clotho.$socket = generateSocketObject();


    function generateSocketObject() {

	    var socket,
		    socketReady,
		    socketQueue = [];

	    var Debugger = new Debug('Socket', '#5555bb');

	    //expecting a string or object
	    function socket_send (data) {
		    if (!socketReady) {
			    Debugger.log('(not ready) queueing request: ', data);
			    socketQueue.push(data);
			    return;
		    }

		    data = angular.isObject(data) ? JSON.stringify(data) : data;

		    Debugger.log('sending data: ', angular.isObject(data) ? data : JSON.parse(data));
		    socket.send(data);
	    }

	    function sendSocketQueue () {
		    Debugger.group('Sending Socket Queue');
		    angular.forEach(socketQueue, function(data, index) {
			    socket_send(data);
		    });
		    Debugger.groupEnd();
		    socketQueue = [];
	    }

	   //usually assume the socket doesn't exist (it would have to be declared manually), but it may if need to connect to a specific port and not what is created here

	    if ($window.$clotho.socket) {
		    socket = $window.$clotho.socket;
	    } else {
		    var pathname = window.location.pathname;

		    var pathDirectory = pathname.substring(0, pathname.lastIndexOf('/'));
		    //check for .html -- if url has hash (e.g. /#!/) then need to strip file name
		    if (/\.html/.test(pathDirectory)) {
			    pathDirectory = pathDirectory.substring(0, pathDirectory.lastIndexOf('/'));
		    }

		    //use protocol wss: for https, default to ws:
		    var socketProtocol = window.location.protocol == 'https:' ? 'wss:' : 'ws:';

		    socket = $window.$clotho.socket = new WebSocket(socketProtocol + "//" + window.location.host + pathDirectory + "/websocket");
	    }

	    if (socket.readyState == 1) {
		    Debugger.log('already exists, sending items in socket Queue...');
		    socketReady = true;
		    sendSocketQueue();
	    } else {
		    socket.onopen = function() {
			    Debugger.log('opened, sending queued items...');
			    socketReady = true;
			    sendSocketQueue();
        }

	    }

	    socket.onerror = function(err) {
		    Debugger.error('socket error', err);
	    };

      socket.onclose = function(evt) {
	      socketReady = false;
        ClientAPI.say({class : "error", text : "Socket Connection Closed", from : "client"});
	      //todo - re-establish connection on loss
      };

	    /*
        //external use - use in other apps, or for custom events
	     //todo - enable as separate service, outside core
        var customChannels = {};

        //external use
        var on = function (eventName, callback) {
            if(!customChannels[eventName]) {
                customChannels[eventName] = [];
            }
            customChannels[eventName].push(callback);
            return [eventName, callback];
        };

        var once = function (eventName, callback) {
            var once_fn = function(args) {
                off([eventName, this]);
                callback(args);
            };
            on(eventName, once_fn);
        };

        var off = function (handle) {
            var t = handle[0];
            customChannels[t] && angular.forEach(customChannels[t], function(idx){
                if(this == handle[1]){
                    console.log("SOCKET\tremoving listener for " + customChannels[t]);
                    customChannels[t].splice(idx, 1);
                }
            });
        };

        //for use internally, in socket.on()
        //future - should this be exposed? Or only reachable by sending message from server?
        var trigger = function (channel, args) {
            customChannels[channel] && angular.forEach(customChannels[channel], function(fn) {
                fn(args);
            });
        };
        */

        /************
         Socket Listener
        ************/

        socket.onmessage = function (obj) {

            obj = JSON.parse(obj.data);

            Debugger.log('received', obj);

            var channel = obj.channel;
            var requestId = obj.requestId;
	          var dataUndefined = angular.isUndefined(obj.data);
            var data = obj.data;

            // it's the ClientAPI method's responsibility to handle data appropriately.
            if (angular.isFunction(ClientAPI[channel])) {
                //Debugger.log("mapping to ClientAPI - " + channel);
                ClientAPI[channel](data);
            }
            /*
            //for custom listeners attached
            else if (typeof customChannels[channel+':'+requestId] == 'function') {
                Debugger.log("mapping to custom listeners - " + channel);
                trigger(channel, data);
            }
            */
            // don't know what to do, so publish to PubSub
            else {
	            //if data field is undefined, send to PubSub.reject to reject the promise.
	            var command = dataUndefined ? 'reject' : 'trigger';
	            if (requestId) {
		            PubSub[command](channel+':'+requestId, data);
	            } else {
		            PubSub[command](channel, data);
	            }
            }
        };

        return {

	          /*
            //For adding custom channels - for use in other apps etc.
            //todo - enable as separate service, outside core
            on : on,
            once : once,
            off : off,
            */

		        state : function () {
			        return socket.readyState;
		        },
            //send a JSON on an arbitrary channel [ repackaged using send() ]
            //note - callback is run on send, not really a callback
            emit: function (channel, data, callback) {
                var packaged = {
                    "channel" : channel,
                    "data" : data
                };
                socket_send(packaged);

                if (typeof callback == 'function') {
	                callback(packaged);
                }
            },

            //send properly packaged and formatted string
            send: socket_send
        }
    }
});
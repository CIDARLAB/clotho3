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

      function createSayError (msg) {
        return {
          class : "error",
          from : "client",
          text : msg
        }
      }

	    //basic check -- most checks should be in communicator
	    function checkMessageValid (msgString) {
		    return msgString.length < 16348;
	    }

	    //expecting a string or object
	    // returns: undefined - put in queue
	    //          true - sent successfully
	    //          false - failed to send
	    function socket_send (data) {
		    //todo - move to socket.readyState
        if (socket.readyState !== 1) {
			    Debugger.log('(not ready) queueing request: ', data);
			    socketQueue.push(data);
			    return;
		    }

		    data = angular.isObject(data) ? JSON.stringify(data) : data;

		    if (checkMessageValid(data)) {
			    Debugger.log('sending data: ', angular.isObject(data) ? data : JSON.parse(data));
			    socket.send(data);
			    return true;
		    } else {
			    Debugger.warn('Message did not pass validation!');
			    ClientAPI.say(createSayError("Message could not be sent to server"));
			    return false;
		    }
	    }

	    function sendSocketQueue () {
		    Debugger.group('Sending Socket Queue');
		    angular.forEach(socketQueue, function(data, index) {
			    socket_send(data);
		    });
		    Debugger.groupEnd();
		    socketQueue = [];
	    }

      function createSocketConnection () {
        var pathname = window.location.pathname;

        var pathDirectory = pathname.substring(0, pathname.lastIndexOf('/'));
        //check for .html -- if url has hash (e.g. /#!/) then need to strip file name
        if (/\.html/.test(pathDirectory)) {
          pathDirectory = pathDirectory.substring(0, pathDirectory.lastIndexOf('/'));
        }

        //use protocol wss: for https, default to ws:
        var socketProtocol = window.location.protocol == 'https:' ? 'wss:' : 'ws:';

        socket = $window.$clotho.socket = new WebSocket(socketProtocol + "//" + window.location.host + pathDirectory + "/websocket");

        socket.onopen = function() {
          Debugger.log('opened, sending queued items...');
          sendSocketQueue();
        };

        socket.onerror = function(err) {
          Debugger.error('socket error', err);
        };

        socket.onclose = function(evt) {
          Debugger.error('socket closed', evt);
          ClientAPI.say(createSayError("Socket Connection Closed (Please Reload)" + (evt.reason ? ' - ' + evt.reason : '') ));

          //todo - actually re-connect with server + ensure events re-attach
          //createSocketConnection();
        };
      }

      //usually assume the socket doesn't exist (it would have to be declared manually), but it may if need to connect to a specific port and not what is created here
      //todo - better check for existence
      if ($window.$clotho.socket) {
        socket = $window.$clotho.socket;
      } else {
        createSocketConnection();
      }

      /************
       Socket Listener
      ************/

      //todo - delegate to service which also handles SSE (once we support them)
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
        //Pubsub if its not in the ClientAPI
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
          state : function () {
            return socket.readyState;
          },

          //send properly packaged and formatted string
          //returns boolean whether message was sent or not
          send: socket_send
      };
    }
});

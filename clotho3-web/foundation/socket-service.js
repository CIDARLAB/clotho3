'use strict';

/* NOTE - see communication dependencies in notes file
    - do socket.on('message') so WebSocket compatible
    - interpret here, send to clothoAPI
*/

Application.Foundation.service('Socket', ['PubSub', 'ClientAPI', function(PubSub, ClientAPI) {

    return (window.$clotho.$socket) ? window.$clotho.$socket : window.$clotho.$socket = generateSocketObject();

    function generateSocketObject() {

        //created in index.html
        var socket = window['$clotho']['socket'];
		
        //internal use (Clotho API) on channel $clotho
        var internal = {};

        //external use - use in other apps, or for custom events
        var customChannels = {};

        /*********
         Socket PubSub
         - Avoid creation of Socket.on() listeners etc. by using custom channels
         *********/

        //INTERNAL USE
        //todo - test these --- may be better to avoid internal registers? just go through PubSub?
        var register = function (eventName, callback) {
            if(!internal[eventName]) {
                internal[eventName] = [];
            }
            internal[eventName].push(callback);
            return [eventName, callback];
        };

        var reg_once = function (eventName, callback) {
            var once_fn = function(args) {
                unregister([eventName, this]);
                callback(args);
            };
            register(eventName, once_fn);
        };

        var unregister = function (handle) {
            var t = handle[0];
            internal[t] && angular.forEach(internal[t], function(idx){
                if(this == handle[1]){
                    console.log("SOCKET\tremoving listener for " + internal[t]);
                    internal[t].splice(idx, 1);
                }
            });
        };

        var internal_trigger = function (channel, args) {
            internal[channel] && angular.forEach(internal[channel], function(fn) {
                fn(args);
            });
        };
        // END INTERNAL USE

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
    /*
     * Fancy socket send with autocallbacks
     */
    socket.idx = -1;
    socket.callbacks = {};
        

    socket.extendedSend = function (message, callback) {
        socket.idx += 1;
        message.requestId = String(socket.idx);
        socket.callbacks[message.requestId] = callback;
        socket.send(JSON.stringify(message));
    };
        /************
         Socket Listener
        ************/

        socket.onmessage = function (obj) {
            
            console.log(obj);
            
            obj = JSON.parse(obj.data);
            var channel = obj.channel;
            var data = obj.data;

        if ("requestId" in obj && obj.requestId in socket.callbacks){
            var callbackKey = obj.requestId;
            var callback = socket.callbacks[callbackKey];
            delete socket.callbacks[callbackKey];
            callback(obj.data)
            return;
        }
            //note - channel reserved for serverAPI
            if (channel == "$clotho") {
                console.log("SOCKET\tchannel $clotho");
                internal_trigger(data.channel, data.data);
                return;
            }

            // it's the ClientAPI method's responsibility to handle data appropriately.
            if (typeof ClientAPI[channel] == 'function') {
                console.log("SOCKET\tmapping to ClientAPI - " + channel);
                ClientAPI[channel](data);
            }
            //for custom listeners attached
            else if (typeof customChannels[channel] == 'function') {
                console.log("SOCKET\tmapping to custom listeners - " + channel);
                trigger(channel, data);
            }
            // don't know what to do, so publish to PubSub
            else {
                console.log("SOCKET\tno listener found for channel: " + channel);
                PubSub.trigger(channel, data);
            }
        };
		   	
		
        return {

            //only for use by Clotho API
            //todo - decide if will be used
            register : register,
            reg_once : reg_once,
            unregister : unregister,


            //For adding custom channels - for use in other apps etc.
            on : on,
            once : once,
            off : off,

            //send a JSON on a 'custom channel' [ repackaged using send() ]
            //future - add callback?
            emit: function (eventName, data) {
                console.log("SOCKET\tdata emitted on channel: " + eventName);

                var packaged = {
                    "channel" : eventName,
                    "data" : data
                };
                socket.send(packaged);
            },
            
            //send properly formatted string on channel message
            send: function (data) {
            	socket.send(data);
        	},

            extendedSend: socket.extendedSend
        }
    }
}]);

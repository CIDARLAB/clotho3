Application.Services.service('PubSubOld', function() {
    /* ABOUT
     // Purpose
        - singleton pub/sub message queue for broadcasting messages between controllers & interacting with server
        - reduce reliance on $rootScope.$broadcast b/c bubbles and slows over large projects
     // PubSub Requirements
        - multiple channels
        - pass in and our arguments (as object?)
        - custom callbacks

     // future - if have a compelling reason to, add listenTo(), stopListening(), listenToOnce() ------- http://backbonejs.org/docs/backbone.html

     //FUTURE - once server communication set up, set up inter-app communication:
     - Notes:
        - some events don't need to be published to the server?
        - Storage events will be handled automatically via localStorage
        - route events should not be published across apps
        - e.g. error saving something wouldn't be meaningful in another page
        - events which should be published:
        - ??
     - Path
        - publish sends messages to server, not to listeners within an app
        - server sends message back to client (so received all pages)
        - relayed into PubSub, all listeners triggered
     - Implement
        - in theory, publish becomes: publish_to_server -> server -> publish_to_client

     // Notes
     - using pub/sub requires calling unsubscribe on $destroy -- need to register watchers to do this

     // Examples
     // More robust - http://www.gridlinked.info/angularJS/modules/MessagingServices.js (jQuery reliance)
     */

    /**
     * future -- decide this later
     * note - combine model_add and model_change
     *
     * API -- STANDARD MESSAGES
     * [messages to be part of a standard API]
     *
     * model_add (uuid, model) -- model added to collector
     * model_destroy (uuid, model) -- model removed from collector
     * model_change (uuid, model) -- model is changed ---- use for listening to any model change
     * model_change:[uuid] (model -- model with <uuid> is changed ---- use for listening for a specific model
     * model_request (uuid) -- model is requested from server
     * model_save (uuid, model) -- model successfully saved to server
     * model_error (uuid, model, xhr) -- saving a model to server fails
     *
     * collector_reset () -- collector is reset (all models deleted)
     *
     * route_change (route, params) -- fired when route changes
     * route:[name] (params) -- fired when a specific route is matched
     *
     */

    //TODO
    // - optimize handling splitting of space-separated events and downstream handling
    // - better handling of model_change AND model_change:[uuid] -- handle splitting or something?

    /***** CONFIG ****/

    //hash of listeners
    var listeners = {};

    //split events passed in by a space
    var eventSplitter = /\s+/;


    /***** HELPERS *****/

    //internal method to add callback to listeners, used in subscribe
    var addToTopic = function(topic, callback) {
        if(!listeners[topic]) {
            listeners[topic] = [];
        }
        listeners[topic].push(callback);
        return [topic, callback];
    };

    //internal method to remove a callback from a given topic
    var removeFromTopic = function(handle) {
        var t = handle[0];

        listeners[t] && angular.forEach(listeners[t], function(idx){
            if(this == handle[1]){
                listeners[t].splice(idx, 1);
            }
        });
    };

    /***** FUNCTIONS ******/

    //testing - log listeners
    var logListeners = function() {
        console.log(listeners);
    };

    /**
     Publish some data on a topic
     @param topic {string} channel to publish on
     @param args {string | object | array}  What the callback for a subscriber would expect.

     note: using function.apply() requires an array be passed, and it's not the fastest. function.call() requires a known number of parameters, and requires a switch depending on the number of parameters.
     note: optimized using invoke(undef) : http://jsperf.com/apply-vs-call-vs-invoke
     */
    var publish = function(topic, args) {
        var topics = topic.split(eventSplitter);
        angular.forEach(topics, function(curTopic) {
            console.log("PUBSUB\tpublish on " + curTopic + " " + args);
            listeners[curTopic] && angular.forEach(listeners[curTopic], function(fn) {
                fn(args);
            });
        });
    };

    /**
     Register a callback to a topic
     @param topic {string} channel to subscribe to. Can pass in multiple, space-separated.
     @param callback {function} handler event. Will be called on publish event to given topic, passed args from publish
     @return handle to pass into unsubscribe. If multiple events are passed in, an array is returned.
     */
    // future - rewrite, handle multiple events better (no if, smarter return)
    var subscribe = function(topic, callback) {
        //handle space-separated events names
        if (eventSplitter.test(topic)) {
            var handles = [];
            var topics = topic.split(eventSplitter);
            angular.forEach(topics, function(curTopic) {
                var handle = addToTopic(curTopic, callback);
                handles.push(handle);
            });
            return handles;
        //for single events passed
        } else {
            return addToTopic(topic, callback);
        }
    };

    /**
     * Register a callback to a topic only a single time
     * @param topic {string}
     * @param callback {function}
     *
     */
    var sub_once = function(topic, callback) {
        var once_fn = function(args) {
            unsubscribe([topic, this]);
            callback(args);
        };
        subscribe(topic, once_fn);
    };

    /**
     Disconnect subscribed function from topic
     @param handle {array} return value from a subscribe() call

     e.g. var handle = PubSub.subscribe( "someTopic", function(){} );
     PubSub.unsubscribe(handle)

     -NOTE: for large PubSub systems, may want to do this another way... right now, loop through array of subscribers on a topic, see if callback matches, splice it out if does. Appears to be how this is commonly done.
     */
    var unsubscribe = function(handle) {
        //handle multiple handles being passed in
        if (angular.isArray(handle[0])) {
            angular.forEach(handle, function(curHandle) {
                removeFromTopic(curHandle);
            });
        }
        else {
            removeFromTopic(handle);
        }
    };

    return {
        logListeners : logListeners,
        publish: publish,
        trigger: publish,
        subscribe: subscribe,
        on : subscribe,
        sub_once: sub_once,
        once : sub_once,
        unsubscribe: unsubscribe,
        off : unsubscribe
    }

});
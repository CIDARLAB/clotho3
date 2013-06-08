Application.Foundation.service('PubSub', function() {

    /* ABOUT
     // Purpose
         - singleton pub/sub message queue for broadcasting messages between controllers & interacting with server
         - reduce reliance on $rootScope.$broadcast b/c bubbles and slows over large projects

     // Notes
        - using pub/sub requires calling unsubscribe on $destroy -- need to register watchers to do this
     */

    /**
     * future -- decide API later
     * note - combine model_add and model_change
     *
     * API -- STANDARD MESSAGES
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

    /*TODO - rewrite with following goals:

     - Track Multiple callbacks for each event
     - Track reference for each callback
     - On() Pass multiple events, space-separated
     - Once() via Flag for Boolean 'once'
     - Each() (unexposed) to handle multiple/single items the same way
     - Appropriate handle return
     - Off() by handle (event, specific callback)
     - Destroy() by reference
     - Clear() by event

     */




    //see if already exists, don't re-instantiate
    return (window.$clotho.$pubsub) ? window.$clotho.$pubsub : window.$clotho.$pubsub = generatePubSubObject();

    function generatePubSubObject() {
        /***** CONFIG ****/

        // hash of listeners
        // listeners[topic] -> [ [fn, flag], [fn, flag], ... ]
        // where flag is true or false... true to delete after run
        var listeners = {};

        //record handles for each reference (i.e. in Angular, $scope.$id
        var references = {};

        //split events passed in by a space
        var eventSplitter = /\s+/;


        /***** HELPERS *****/

        //track listeners for a given reference
        var addToRef = function(ref, handle) {
            if(!references[ref]) {
                references[ref] = [];
            }
            references[ref].push(handle);
        };

        //internal method to add callback to listeners, used in subscribe
        var addToTopic = function(topic, callback, ref, one) {
            if(!listeners[topic]) {
                listeners[topic] = [];
            }

            listeners[topic].push([callback, one]);

            var handle = [topic, callback];
            if (ref && typeof ref != 'undefined') {
                addToRef(ref, handle);
            }
            return handle;
        };

        //internal method to remove a callback from a given topic
        var removeFromTopic = function(handle) {

            var t = handle[0];

            listeners[t] && angular.forEach(listeners[t], function(array, idx){
                if(array[0] == handle[1]){
                    //console.log("splicing - " + t);
                    listeners[t].splice(idx, 1);
                } else {
                    //console.log("not splicing: " + t);
                    //console.log(array[0]);
                    //console.log(handle[1]);

                }
            });
        };

        //testing - log listeners
        var logListeners = function() {
            console.log(listeners);
        };


        /***** FUNCTIONS ******/

        /**
         @name PubSub.trigger
         @description
         Publish some data on a topic
         @param topic {string} channel to publish on, can be multiple space-separated
         @param args {string | object | array}  What the callback for a subscriber would expect.

         note: using function.apply() requires an array be passed, and it's not the fastest. function.call() requires a known number of parameters, and requires a switch depending on the number of parameters.
         note: optimized using invoke(undef) : http://jsperf.com/apply-vs-call-vs-invoke
         */
        var trigger = function(topic, args) {
            var topics = topic.split(eventSplitter);
            angular.forEach(topics, function(curTopic) {
                // testing
                // console.log("PUBSUB\tpublish on " + curTopic + " " + args);
                listeners[curTopic] && angular.forEach(listeners[curTopic], function(array, idx) {
                    array[0](args);
                    if (array[1]) {
                        //console.log("splicing inside trigger: " + curTopic);
                        //console.log(array[0]);
                        listeners[curTopic].splice(idx, 1);
                    } else {
                        //console.log("not splicing in trigger: " + curTopic);
                    }
                });
            });
        };

        /**
         @name PubSub.on
         @description
         Register a callback to a topic

         @param topic {string} channel to subscribe to. Can pass in multiple, space-separated.
         @param callback {function} handler event. Will be called on publish event to given topic, passed args from publish
         @param ref {string} $scope.$id (or other name), which when destroyed should delete its associated listeners. Avoid passing in objects.
         @param one {boolean} Flag to run the callback only once
         @return handle to pass into unsubscribe. If multiple events are passed in, an array is returned.
         */
        // future - rewrite, handle return multiple events better (no if, smarter return)
        var on = function(topic, callback, ref, one) {
            ref = ref || null;
            one = one || false;

            //handle space-separated events names
            if (eventSplitter.test(topic)) {
                var handles = [];
                var topics = topic.split(eventSplitter);
                angular.forEach(topics, function(curTopic) {
                    var handle = addToTopic(curTopic, callback, ref, one);
                    handles.push(handle);
                });
                return handles;
            } else {
                //console.log(callback);
                return addToTopic(topic, callback, ref, one);
            }
        };

        /**
         * @name PubSub.once
         * @description
         * Register a callback to a topic only a single time
         * @param topic {string}
         * @param callback {function}
         * @param ref {string} $scope.$id (or other name), which when destroyed should delete its associated listeners. Avoid passing in objects.
         *
         */
        var once = function(topic, callback, ref) {
            on(topic, callback, ref, true);
        };

        /**
         @name PubSub.off
         @description
         Disconnect subscribed function from topic
         @param handle {array} return value from a subscribe() call

         e.g. var handle = PubSub.subscribe( "someTopic", function(){} );
         PubSub.unsubscribe(handle)

         -NOTE: for large PubSub systems, may want to do this another way... right now, loop through array of subscribers on a topic, see if callback matches, splice it out if does. Appears to be how this is commonly done.
         */
        var off = function(handle) {
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

        /**
         * @name PubSub.destroy
         *
         * @param ref {string|object} $scope.$id (or other string) which when destroyed should remove its associated listeners. Can pass object, but slower so should be avoided.
         *
         * @description
         * Removes all listeners for an associated reference (i.e. in Angular, $scope.$id)
         * E.g. $scope.$on('$destroy', Clotho.silence($scope.$id));
         */
        var destroy = function(ref) {

            //console.log(references);
            //console.log(ref);

            var ref_handles = references[ref];
            //console.log(ref_handles);

            angular.forEach(ref_handles, function(handle) {
                //console.log(handle);
                removeFromTopic(handle);
            });

            //clean up ref map
            references[ref] = [];
        };

        /**
         * @name PubSub.clear
         *
         * @param topic {string} Topic to clear, can be space-separated list
         *
         * @description
         * Clears all listeners a topic.
         * Clears all listeners if "all" is passed.
         * Does nothing if no topic passed.
         */
        var clear = function(topic) {
            if (topic == "all") {
                listeners = {};
                references = {};
            }

            var topics = topic.split(eventSplitter);
            angular.forEach(topics, function(curTopic) {
                listeners[curTopic] = [];
            });

            //todo - clean up references array also --- this is tricky

        };


        return {
            logListeners : logListeners,
            trigger: trigger,
            on : on,
            once : once,
            off : off,
            destroy : destroy,
            clear : clear
        }
    }
});
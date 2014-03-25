/**
 * @name clothoAPI
 *
 * @description
 * This is the server Clotho API - commands that will be issued to the server
 *
 * Note - General Philosophy and Usage
 *  - These methods are exposed to users for their use. They should be easy to use.
 *  - These methods are meant to communicate with the socket
 *  - They will be exposed outside of Angular
 *
 *  Most functions are written so that resulting behavior is not captured in a return, but rather by behavior within the ClientAPI. For some cases, returns are necessary for the function to make sense (e.g. Clotho.get() ).
 *
 *  At present, the goal is write them such that all resulting behavior will be captured by the messages returned over the socket from the server -- i.e. behaviour within clientAPI, PubSub, etc. This may not be possible, but regardless, the channel '$clotho' is reserved for the serverAPI (and falls through on clientAPI).
 * The main problem is get()... which by its nature should return the object to the callback. However, it may make sense to encapsulate this functionality in the collector. See note below
 *
 * Future todos
 *  - return angular.noop for empty callbacks, or deal with better
 */

angular.module('clotho.core').service('Clotho',
	function(Socket, Collector, PubSub, Debug, $window, $q, $rootScope, $location, $timeout) {

//attach to window, only instantiate once
return ($window.$clotho.api) ? $window.$clotho.api : $window.$clotho.api = generateClothoAPI();

function generateClothoAPI() {

    /**********
       Set up
     **********/
    var Debugger = new Debug('ClothoAPI', '#cc5555');

    //socket communication
    var fn = {};

		//note that angular.toJson will strip $-prefixed keys, so should be avoided
    fn.send = function(pkg) {
        Socket.send(angular.toJson(pkg));
    };
    fn.pack = function(channel, data, requestId, options) {
        return {
            "channel" : channel,
            "data" : data,
            "requestId" : requestId,
            "options" : options
        };
    };
    fn.emit = function(eventName, data, requestId, options) {
        fn.send(fn.pack(eventName, data, requestId, options));
    };

    //helper functions

		//pending ES6

		var numberAPICalls = new (function () {
			var commandNum = 0;
			this.next = function () {
				commandNum += 1;
				return commandNum.toString();
			};
		});

	/**
	 * @name fn.emitSubCallback
	 * @param channel {String} API Channel
	 * @param data {*} API Data
	 * @param func {Function} Callback Function, run when data received from server
	 * @param options {Object} API Options
	 * @returns {Promise} Resolved on receipt of data from server, rejected after 5 second timeout
	 */
	fn.emitSubCallback = function (channel, data, func, options) {
		var deferred = $q.defer(),
			requestId = Date.now().toString() + numberAPICalls.next();

		if (!angular.isFunction(func)) {
			func = angular.noop;
		}

		//timeout our requests at 5 seconds
		var timeoutPromise = $timeout(function () {
			deferred.reject(null)
		}, 5000);

		//todo - reject promise if no data is sent over socket
		PubSub.once(channel + ':' + requestId, function (data) {
			$timeout.cancel(timeoutPromise);
			deferred.resolve(data);
			func(data);
		}, '$clotho');

		fn.emit(channel, data, requestId, options);
		return deferred.promise;
	};

	//if simply want a response and resolve promise, no further logic required
	fn.emitSubOnce = function (channel, data, options) {
		return fn.emitSubCallback(channel, data, angular.noop, options);
	};



    /**********
       Config
     **********/


    /**********
      Testing
     **********/


    /**********
         API
     **********/

    /**
     * @name Clotho.login
     *
     * @param username
     * @param password
     *
     * @description
     * Authenticate with Clotho
     *
     * @returns {Promise} result of login
     */
    var login = function clothoAPI_login(username, password) {
        var cred = {username: username, password: password};
        return fn.emitSubOnce('login', cred);
    };

		/**
		 * @name Clotho.logout
		 *
		 * @description
		 * Logout of Clotho
		 *
		 * @returns {Promise} result of login
		 */
		var logout = function clothoAPI_logout() {
			return fn.emitSubOnce('logout', '');
		};

    /**
     * @name Clotho.get
     *
     * @param {string|array} uuid UUID of Sharable, or an array of string UUIDs
     * @param {object} options
     * @param {boolean=} synchronous   Default false for async, true to return value synchronously
     *
     * @description
     * Request a Sharable from the server. Returns the object description for the selected objects. If a selector is ambiguous, clotho will return the first one returned by MongoDB - which is basically arbitrary. Clotho will send an error message on the 'say' channel if an object could not be retrieved.
     *
     * @returns {Promise | object} Promise to be fulfilled with JSON representation of Sharable, or the object for synchronous (blocking).
     * NB: when a promise is returned, it is resolved on $$v if you don't use .then()....
     *
     * @usage
     * RECOMMENDED
     * Clotho.get(uuid).then(function (result) {
     *  var workWithMe = result;
     * });
     *
     * ALTERNATIVELY
     * var workWithMe = Clotho.get(uuid);
     * ... later ...
     * var populated = workWithMe.$$v
     *
     */
    var get = function clothoAPI_get(uuid, options, synchronous) {

        //default to false (return promise)
        synchronous = typeof synchronous != 'undefined' ? !!synchronous : false;

        return synchronous ? get_sync(uuid, options) : get_async(uuid, options);
    };

    var get_async = function clothoAPI_get_async(uuid, options) {

		    if (angular.isUndefined(uuid)) {
			    return $q.when();
		    }

	      //check collector
	      var retrieved = Collector.retrieveModel(uuid);

	      if (!!retrieved) {
	          var deferred = $q.defer();
	          deferred.resolve(retrieved);
	          return deferred.promise;
	      } else {
	          var callback = function(data) {
	              Collector.storeModel(uuid, data);
	          };

	          return fn.emitSubCallback('get', uuid, callback, options);
	      }
    };

    var get_sync = function clothoAPI_get_sync(uuid, options) {

		    if (angular.isUndefined(uuid)) {
			    return;
		    }

        //check collector
        var data = Collector.retrieveModel(uuid);

        if (!!data) {
            return data
        } else {
            //future - need REST API

            /* OLD ATTEMPT
             var value = {};
             value.$then = deferred.promise.then(function (response) {
                var then = value.$then;
                if (response) {
                    angular.copy(response, value);
                    value.$then = then;
                }
             });
             return value;
             */
        }
    };

    /**
     * todo - update doc
     *
     * @name Clotho.set
     *
     * @param {object} sharable  JSON representation of Sharable containing new data - may be partial (just a few fields), but must contain an ID or UUID
     *
     * @description
     * Updates a sharable with passed object. Sets the fields present in the spec to the values in the spec. Fields missing from the spec are unchanged. If a field is missing from the specified object, it will be created. Clotho will send an error message on the 'say' channel if an update could not be applied.
     *
     * Server will emit a collect(uuid) call, upon object being updated
     */
    var set = function clothoAPI_set(sharable) {

        if (angular.isEmpty(sharable)) return false;

        //strip $$v (angular promise value wrapper)
        while (sharable.$$v) {
            sharable = sharable.$$v;
        }

        var callback = function() {
            Collector.storeModel(sharable.id, sharable);
        };

        return fn.emitSubCallback('set', sharable, callback);
    };


    /**
     * @name Clotho.clone
     *
     * @param {string} uuid
     *
     * @description
     * copy a sharable, which won't be updated later
     *
     * @return
     * returns a full JSON representation of a sharable
     */

    var clone = function clientAPIClone(uuid) {

    };

	/**
	 * @name Clotho.watch
	 *
	 * @param {string} uuid UUID of Sharable to watch for changes
	 * @param {object|function} action if Object, object to be updated (using angular.extend) using passed model for given UUID. if Function, function to run on change, passed the new value, with this equal to the reference passed
	 * @param {string} reference Reference (e.g. $scope) for element to unlink listener on destroy. Passing in a $scope object (e.g. from a controller or directive) will automatically handle deregistering listeners on its destruction.
	 * @param {boolean} overwriteExistingObj If truthy, will extend an empty object, removing existing fields. Only applies if extending an object
	 *
	 * @description
	 * Watches for published changes for a given uuid, updating the object using angular.extend
	 */
		var watch = function clothoAPI_watch(uuid, action, reference, overwriteExistingObj) {
			reference = typeof reference != 'undefined' ? reference : null;
			PubSub.on('update:'+uuid, function(model) {
				if (angular.isFunction(action)) {
					action.apply(reference, model);
				}
				else {
					if (!!overwriteExistingObj) {
						angular.extend({}, model);
					} else {
						angular.extend(action, model);
					}
				}
			}, reference);
		};


    /**
     * @name Clotho.listen
     *
     * @param {string} channel PubSub Channel to listen to
     * @param {function} callback Function to be executed on change
     * @param {string=|object=} reference Reference (e.g. $scope) for element to unlink listener on destroy. Passing in a $scope object (e.g. from a controller or directive) will automatically handle deregistering listeners on its destruction.
     *
     * @description
     * Watches for published events on a given channel, and runs a callback on the event
     */
    var listen = function clothoAPI_listen(channel, callback, reference) {
        reference = typeof reference != 'undefined' ? reference : null;

        PubSub.on(channel, function clothoAPI_listen_callback(data) {
            callback(data);
        }, reference);
    };

    /**
     * @name Clotho.trigger
     * @param channel {string}
     * @param data {*=}
     * @description
     * Publish an event directly to PubSub, bypassing the socket and server.
     */
    var trigger = function clothoAPI_trigger(channel, data) {
        PubSub.trigger(channel, data);
    };


    /**
     * @name Clotho.silence
     *
     * @param {string} reference Reference used when registering listeners. Will destroy all listeners associated with this reference.
     * Passing "all" destroys all listeners.
     *
     * @description
     * Destroys listener functions associated with a given reference
     *
     */
	   //todo - deprecate
    var silence = function clothoAPI_silence(reference) {
        PubSub.destroy(reference);
    };

    /**
     * @name Clotho.emit
     *
     * @param {string} channel
     * @param {object} args Sends {} if null
     *
     * @description
     * Emit an object on a custom channel message to the server
     */
    var emit = function clothoAPI_emit(channel, args) {
        fn.emit(channel, args || {});
    };

    /**
     * @name Clotho.broadcast
     *
     * @param {string} channel
     * @param {object} args
     *
     * @description
     * Broadcast an object on a given channel on the client
     */
    var broadcast = function clothoAPI_broadcast(channel, args) {
        fn.emit('broadcast', fn.pack(channel, args));
    };

    /**
     * @name Clotho.on
     * @param channel
     * @param callback
     * @returns {array} handle
     *
     * @description
     * Registers a listener directly with the socket, bypassing PubSub
     */
    var on = function clothoAPI_on(channel, callback) {
        return Socket.on(channel, callback);
    };

    /**
     * @name Clotho.once
     * @param channel
     * @param callback
     * @returns {array} handle
     *
     * @description
     * Registers a listener ONCE directly with the socket, bypassing PubSub
     */
    var once = function clothoAPI_once(channel, callback) {
        return Socket.once(channel, callback);
    };

    /**
     * @name Clotho.off
     * @param handle Handle returned by Clotho.on() or Clotho.once()
     *
     * @description
     * Removes a listener that had been attached to the Socket
     */
    var off = function clothoAPI_off(handle) {
        Socket.off(handle);
    };



    /**
     * @name Clotho.query
     *
     * @param obj Constraints for query.
     * @param options Options, e.g. fields to return
     *
     * @example To get all schemas, Clotho.query({"schema" : "Schema"})
     *
     * @description
     * Returns all objects that match the fields provided in the spec.
     */
    var query = function(obj, options) {
        var callback = function(data) {
            //store query
            //todo

            //store models
            //future - when not sending whole model, extend what exists
            angular.forEach(data, function(sharable) {
                Collector.storeModel(sharable.id, sharable);
            })
        };
        return fn.emitSubCallback('query', obj, callback, options);
    };


    /**
     * @name Clotho.create
     *
     * @param {object} obj A json object or list of json objects that describe(s) an instance or instances of a clotho schema.
     *
     * @description
     * Creates the described objects. Returns a list of uuids of the created objects, in the same order as they were sent in. If an object was not created, that object's entry in the list of uuids will be null. Reasons for failing to create the failed objects will be provided in the 'say' channel.
     *
     * @returns {object} The object if created successfully, null otherwise
     */
    var create = function clothoAPI_create(obj) {
        return fn.emitSubOnce('create', obj);
    };

    /**
     *
     * todo - update doc
     *
     * @name Clotho.destroy
     *
     * @param {string|array} uuid UUID of Sharable, or an array of string UUIDs
     *
     * @description
     * Destroys the listed objects. Does not destroy anything if the selector is ambiguous. Clotho will send an error message on the 'say' channel if destroy fails for an object.
     *
     */
    var destroy = function clothoAPI_destroy(uuid) {

        var callback = function() {
            Collector.removeModel(uuid);
        };

        return fn.emitSubCallback('destroy', uuid, callback);
    };

    /**
     * @name Clotho.edit
     *
     * @param {string} uuid
     *
     * @description
     * Open a sharable in an editor
     *
     */
    var edit = function clothoAPI_edit(uuid) {
        $location.path("/editor/" + uuid);
    };

    /**
     * @name Clotho.revert
     *
     * @param {string} uuid UUID of resource with unwanted changes
     * @param {string} timestamp Date of desired version
     *
     * @description
     * Revert a model to an earlier version
     *
     * @returns {object} version at timestamp of resource with passed UUID
     * todo - force collection of new model
     *
     */
    var revert = function clothoAPI_revert(uuid, timestamp) {
        var packaged = {"uuid" : uuid, "timestamp" : timestamp};
        return fn.emitSubOnce('revert', packaged);
    };


    /**
     * @name Clotho.validate
     *
     * @param obj
     *
     * @description
     * Passes an object or array of objects to the server to be validated, according to the rules of the object(s)'s respective schema(s)
     *
     * @returns {array} array of results of validation. Empty object means object passed validation. Populated object encountered problems. Fields are listed which were problematic, with error messages.
     */
    var validate = function (obj) {
        return fn.emitSubOnce('validate', obj);
    };

    /**
     * @name Clotho.show_old
     *
     * @param {string} uuid
     * @param {object=} args
     *  - ? [custom]
     *  - position object
     *      - parent : { uuid : [uuid of div to insert into] }
     *      - absposition : { {int} x, {int} y }
     *
     * @description
     * Request a GUI from the server. JSON will be validated and pushed to the client, and necessary resources will be collected via a series of client.collect() calls, then displayed via client.show().
     *
     * @notes
     * viewID should be kept on the server, not passed explicitly to function
     *
     */
    var show_old = function clothoAPI_show(uuid, args) {
        var packaged = {
            "uuid" : uuid,
            "args" : args
        };
        fn.api.emit('show_old', packaged);
    };


    /**
     * @name Clotho.show
     *
     * @param {string} viewId
     * @param {object} options Example given below
     *
     {
         "target" : <css selector>   // suggested, otherwise placed outside ng-view
     }
     */
    var show = function(viewId, options) {
	    var packaged = {
		    "viewId" : viewId,
		    "options" : options
	    };
        fn.emit('show', packaged);
    };


    /**
     * @name Clotho.share
     *
     * @description Opens a modal to share the current page
     */
    var share = function() {
	      //todo - need to set up share (outside API)
	      console.log('need to set up share');
        //$modal.share($location.absUrl());
    };


    /**
     * @name Clotho.log
     *
     * @param {string} msg
     *
     * @description
     * Writes a message to the Server log for the Doo managing the current process
     *
     */
    var log = function clothoAPI_log(msg) {
        fn.emit('log', msg);
    };

    /**
     * @name Clotho.say
     *
     * @param {string} msg
     * @param {string} userID
     *
     * @description
     * Post a message to the Command Bar dialog, or to another user's dialog
     *
     */
    var say = function clothoAPI_say(msg, userID) {
        var packaged = {
            "userID" : userID || "",
            "msg" : msg,
            "timestamp" : Date.now()
        };
        fn.emit('say', packaged);
    };

    /**
     * @name Clotho.notify
     *
     * @param {object} data
     *
     * @description
     * notify the server of an event (e.g. model change in collector)
     *
     */
    var notify = function clothoAPI_notify(data) {
        fn.emit('notify', data);
    };

    /**
     * @name Clotho.alert
     *
     * @param {string} msg
     * @param {string} userID
     *
     * @description
     * Post an alert to the current user, or to another user
     *
     */
    var alert = function clothoAPI_alert(msg, userID) {
        var packaged = {
            "userID" : userID,
            "msg" : msg
        };
        fn.emit('alert', packaged);
    };

    var autocomplete = function(query,  options) {
      //todo - use $cacheFactory to cache searches

	    //catch all, publish on autocomplete for now
	    var callback = function(data) {
		    //todo - only publish if newer
		    PubSub.trigger('autocomplete', data)
	    };
	    var packaged = {
        "query" : query
      };
      return fn.emitSubCallback('autocomplete', packaged, callback, options);
    };

    var autocompleteDetail = function(uuid) {

        var deferred = $q.defer();
        //check the collector first
        if (Collector.hasItem('detail_'+uuid)) {
            deferred.resolve(Collector.retrieveModel('detail_'+uuid));
        } else {
            var packaged = {
                "uuid" : uuid
            };
            fn.emit('autocompleteDetail', packaged);

            //testing
            //PubSub.once('autocompleteDetail_'+'function_id123', function(data) {
            PubSub.once('update:detail_'+uuid, function(data) {
                deferred.resolve(data);
            }, '$clotho');
        }

        return deferred.promise;
    };





    var submit = function(query) {
        return fn.emitSubOnce('submit', query);
    };

    /**
     * @name Clotho.run
     *
     * @param {string} func Object UUID or name (as name or module.name) indicating the function to be run (follows get semantics for ambiguous selectors)
     * @param {object} args A JSON object with key-value pairs providing the argument values to the function.
     *
     * @description
     * Runs a function on the server on supplied arguments.
     * Clotho will send an error message on the 'say' channel if there is an error during function execution, there is no function matching the specifier, or there are ambiguously specified arguments.
     */
    var run = function clothoAPI_run(func, args) {
        var packaged = {
            "id" : func,
            "args" : args
        };

        return fn.emitSubOnce('run', packaged);
    };


    /**
     * @name Clotho.recent_deprecated
     * @note This is the old implementation. Sends request on requestRecent, and expects response on displayRecent
     * @returns {Promise}
     */
    var recent_deprecated = function() {
        fn.emit('requestRecent', {});

        var deferred = $q.defer();

        PubSub.once('displayRecent', function(data) {
            deferred.resolve(data);
        }, 'clothoAPI');

        return deferred.promise;
    };

    /**
     * @name Clotho.recent
     *
     * @param {Object=} params Restrictions to your query
     *
     * //todo - examples
     *
     * @description
     * Request your most recently / commonly used sharables
     *
     * @returns {Promise} Data once returned from server
     */
    var recent = function(params) {
        return fn.emitSubOnce('recent', params);
    };


    // ---- TO BE IMPLEMENTED LATER ----

    /**
     * @name Clotho.learn
     *
     * @param {string} query
     * @param {string} assoc
     *
     * @description
     * Learn to associate a string (query) with another string / command (assoc)
     */
    var learn = function clothoAPI_learn(query, assoc) {
        var packaged = {
            "query" : query,
            "association" : assoc
        };
        return fn.emitSubOnce('learn', packaged);
    };

    /**
     * @name Clotho.startTrail
     *
     * @param {string} uuid
     *
     * @description
     * start a trail with a given uuid
     */
    var startTrail = function clothoAPI_startTrail(uuid) {
        $location.path("/trails/" + uuid);
    };

    /**
     * @name Clotho.gradeQuiz
     *
     * @param {*} questionValue
     * @param {*} input
     * @param {string} answerGen ID of function to run to generate answer
     *
     * @description
     * wrapper for grade quiz funciton on server (easier to change in one place)
     *
     */
    var gradeQuiz = function clothoAPI_gradeQuiz(questionValue, input, answerGen) {
        return run('gradeQuiz', [questionValue, input, answerGen]);
    };


    return {
        //api
        login : login,
        logout : logout,
        get : get,
        set : set,
        query : query,
        create : create,
        edit : edit,
        validate : validate,
        revert : revert,
        destroy : destroy,
        show : show,
        show_old : show_old,
        say : say,
        log : log,
        alert : alert,
        run : run,
        recent: recent,
        notify : notify,
        gradeQuiz : gradeQuiz,

		    //searchbar
		    submit: submit,
		    autocomplete : autocomplete,
		    autocompleteDetail : autocompleteDetail,

        //toolkit
        watch : watch,
        listen : listen,
        silence : silence,
        trigger: trigger,
        emit : emit,
        broadcast : broadcast,
        on : on,
        once : once,
        off : off,
        share : share

    }

}
});

'use strict';

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
 *  Most functions are written so that resulting behavior is not captured in a return, but rather by behavior within the ClientAPI. For some cases, returns are necessary for the function to make sense (e.g. Clotho.get()). Most
 *
 *  At present, the goal is write them such that all resulting behavior will be captured by the messages returned over the socket from the server -- i.e. behaviour within clientAPI, PubSub, etc. This may not be possible, but regardless, the channel '$clotho' is reserved for the serverAPI (and falls through on clientAPI).
 * The main problem is get()... which by its nature should return the object to the callback. However, it may make sense to encapsulate this functionality in the collector. See note below
 *
 * Note - Dependencies & the Socket
 * This service is kept inside of Angular to take advantage of Angular core services. It should have no dependencies on services etc. written by us, in order to avoid circular dependencies.
 * While some of the functionality encapsulated by the Socket service would certaintly be useful here, and more importantly the socket.on() listener must be duplicated, that is of less concern than the serverAPI being exposed to other services (esp. Collector) and the functionality gained by keeping it within Angular.
 *
 * Note - Collector usage
 * In the app, retrieval of models should be handled by the Collector. The function can return a promise and fulfill it, but dependent on two different services. It could be handled as two steps ( in retrieveModel() ):
 *      (1) Clotho.get(uuid) to get the object, then
 *      (2) PubSub.sub_once(model_change, uuid) to set the object
 *
 * Note - Alternative Setups
 * (1) ServerAPI dependent on Socket. All other services interact only with non-Angular version.
 *      - If hacking this to avoid Angular DI, that's probably not a good idea.
 *      - by adding dependency, no other service can list this one as a dependency
 *      - Doesn't functionally accomplish much, except cleaner code
 *      + serverAPI could also interface with PubSub
 *      + only one instance of socket.on() listener
 *
 * Future todos
 *  - return angular.noop for empty callbacks, or deal with better
 */

Application.Foundation.service('Clotho', ['Socket', 'Collector', 'PubSub', '$q', '$rootScope', '$location', '$timeout', function(Socket, Collector, PubSub, $q, $rootScope, $location, $timeout) {

    /**********
       Set up
     **********/
    var fn = {};

    fn.send = function(pkg) {
        Socket.send(JSON.stringify(pkg));
    };
    fn.pack = function(channel, data) {
        return {
            "channel" : channel,
            "data" : data
        };
    };
    fn.emit = function(eventName, data) {
        fn.send(fn.pack(eventName, data));
    };

    fn.api = {};
    fn.searchbar = {};

    fn.api.emit = function (eventName, data) {
        fn.emit("api",
            fn.pack(eventName, data)
        )
    };
    fn.searchbar.emit = function (eventName, data) {
        fn.emit("searchbar",
            fn.pack(eventName, data)
        )
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
     * @name Clotho.get
     *
     * @param {string} uuid UUID of Sharable
     * @param {boolean=} synchronous   Default false for async, true to return value synchronously
     *
     * @description
     * Request a Sharable from the server
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
     * todo - add callback? Promises too complicated?
     * todo - look at angular resource - doesn't populate on $$V???
     */
    var get = function clothoAPI_get(uuid, synchronous) {

        //default to false (return promise)
        synchronous = typeof synchronous != 'undefined' ? synchronous : false;

        var deferred = $q.defer();

        //check collector
        var data = Collector.retrieveModel(uuid);

        if (data && typeof data != 'undefined') {
            if (synchronous) {
                return data;
            } else {
                deferred.resolve(data);
                //deferred.promise.then(function(result) {deferred.promise = result});
                return deferred.promise;
            }
        }
        //otherwise request over socket, will be collected, listen for event
        else {
            fn.api.emit('get', uuid);

            PubSub.once('model_change:'+uuid, function(data) {
                //console.log("CLOTHOAPI\t");
                //console.log(data);
                $rootScope.$safeApply(deferred.resolve(data));
            }, '$clotho');

            if (!synchronous) {
                return deferred.promise;
            } else {
                //fixme - need to timeout or something for synchronous use -- look for $timeout.resolve solution

                var value = {};
                value.$then = deferred.promise.then(function (response) {
                    var then = value.$then;
                    if (response) {
                        angular.copy(response, value);
                        value.$then = then;
                    }
                }).then;
                return value;
            }
        }
    };

    //for testing
    var get_url = function clothoAPI_get_url(uuid) {
        var deferred = $q.defer();

        var data = Collector.retrieveModel(uuid);
        if (data && typeof data != 'undefined') {
            deferred.resolve(data);
            return deferred.promise;
        } else {
            fn.api.emit('get_url', uuid);

            PubSub.once('model_change:'+uuid, function(data) {
                console.log("got template url for uuid " + uuid);
                $rootScope.$safeApply(deferred.resolve(data));
            });

            return deferred.promise;
        }
    };

    //for testing
    var get_template = function clothoAPI_get_template(uuid) {
        fn.api.emit('get_template', uuid);
    };

    //for testing
    var get_script = function clothoAPI_get_script(uuid) {
        fn.api.emit('get_script', uuid);
    };

    /**
     * @name Clotho.set
     *
     * @param {string} uuid UUID of Sharable to alter
     * @param {object} data  JSON representation of new data, may be partial (just a few fields)
     *
     * @description
     * Updates a sharable with passed object. The passed object may only contain a few fields of the Sharable, but must pass validation on the server.
     *
     * Server will emit a collect(uuid) call, upon object being updated
     */
    var set = function clothoAPI_set(uuid, data) {
        //strip $$v in case use promise to set multiple times
        while (data.$$v) {
            data = data.$$v;
        }

        var packaged = {
            "uuid" : uuid,
            "data" : data
        };
        fn.api.emit('set', packaged);
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

    //if pass in scope, handle removing listeners
    function silenceRefIfScope(reference) {
        if (reference.$evalAsync && reference.$watch) {
            reference.$on('$destroy', function() {
                //console.log("silencing " + reference.$id);
                silence(reference);
            })
        }
    }

    /**
     * @name Clotho.watch
     *
     * @param {string} uuid UUID of Sharable to watch for changes
     * @param {function} callback Function to be executed on change
     * @param {string=} reference Reference (e.g. $scope.$id) for element to unlink listener on destroy. Should be included, or will bloat.
     *
     * @description
     * Watches for published changes for a given uuid
     *
     * todo - smarter reference logic
     */
    var watch = function clothoAPI_watch(uuid, callback, reference) {
        reference = typeof reference != 'undefined' ? reference : null;
        PubSub.on('model_change:'+uuid, function(model) {
            $rootScope.$safeApply(callback(model));
        }, reference);

        silenceRefIfScope(reference);
    };


    /**
    // todo - mix into watch, check for function || array
     * @name Clotho.watch2
     *
     * @param {string} uuid UUID of Sharable to watch for changes
     * @param {object} scope Object to be updated
     * @param {string} field field in object to be updated
     * @param {string} reference Reference (e.g. $scope.$id) for element to unlink listener on destroy
     *
     * @description
     * Watches for published changes for a given uuid
     * note: the object and field cannot be passed in together due to the way javascript handles references. Passing them in separately allows a given field to be changed in the object itself, and the change to exist beyond the scope of this function.
     *
     */
    var watch2 = function clothoAPI_watch2(uuid, scope, field, reference) {
        reference = typeof reference != 'undefined' ? reference : null;
        PubSub.on('model_change:'+uuid, function clothoAPI_watch2_callback(model) {
            $rootScope.$safeApply(scope[field] = model);
        }, reference);

        silenceRefIfScope(reference);
    };

    /**
     * @name Clotho.listen
     *
     * @param {string} channel PubSub Channel to listen to
     * @param {function} callback Function to be executed on change
     * @param {string=} reference Reference (e.g. $scope) for element to unlink listener on destroy. Passing in a $scope object (e.g. from a controller or directive) will automatically handle deregistering listeners on its destruction.
     *
     * @description
     * Watches for published events on a given channel, and runs a callback on the event
     */
    var listen = function clothoAPI_listen(channel, callback, reference) {
        reference = typeof reference != 'undefined' ? reference : null;

        PubSub.on(channel, function clothoAPI_listen_callback(data) {
            $rootScope.$safeApply(callback(data));
        }, reference);

        silenceRefIfScope(reference);
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
    var silence = function clothoAPI_silence(reference) {
        PubSub.destroy(reference);
    };

    /**
     * @name Clotho.emit
     *
     * @param {string} channel
     * @param {object} args
     *
     * @description
     * Emit an object on a custom channel message to the server
     */
    var emit = function clothoAPI_emit(channel, args) {
        fn.api.emit(channel, args);
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
        fn.api.emit('broadcast', fn.pack(channel, args));
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




    var query = function(obj) {
        fn.searchbar.emit("query", obj);
    };

    /**
     * @name Clotho.create
     *
     * @param {string | object} schema
     * - string: UUID of Schema to use for validation
     * - object: JSON representation of Schema to be used
     *
     * @param {object} data Model for the resource to be created
     * @param {function()=} callback
     *
     * @description
     * Create a Sharable on the server
     *
     * @returns {object} The object if created successfully, null otherwise
     */
    var create = function clothoAPI_create(schema, data, callback) {
        var packaged = {
            "schema" : schema,
            "data" : data
        };
        fn.api.emit('create', packaged);
        //todo - return
    };

    /**
     * @name Clotho.destroy
     *
     * @param {string} uuid
     *
     * @description
     * Delete a sharable, provided you have the given permissions
     *
     */
    var destroy = function clothoAPI_destroy(uuid) {
        fn.api.emit('destroy', uuid);
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
        //fn.api.emit('edit', uuid);

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
     *
     */
        //todo - callback?
    var revert = function clothoAPI_revert(uuid, timestamp) {
        fn.api.emit('revert', {"uuid" : uuid, "timestamp" : timestamp});
    };

    /**
     * @name Clotho.show
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
    var show = function clothoAPI_show(uuid, args) {
        var packaged = {
            "uuid" : uuid,
            "args" : args
        };
        fn.api.emit('show', packaged);
    };


    //testing
    var show_simple = function(obj) {
        fn.api.emit('show_simple', obj);
    };

    /**
     * @name Clotho.bootstrap
     *
     * @param appInfo {object} Object with necessary information to bootstrap, minimally including:
     * {
     *      "moduleName" : <Name of module as defined in Angular>
     *      "moduleUrl" : "<URL to module js file>",
     * }
     *
     * @returns {array} Selectors in form: [<appUUID>, <jQuery Selector>]
     * @description
     * Load a widget and bootstrap it
     * note - jQuery dependency
     */
    var nextWidgetUUID = 0;
    var bootstrap = function clothoAPI_bootstrap(appInfo) {

        nextWidgetUUID++;
        var appUUID = nextWidgetUUID;

        var deferred = $q.defer();

        /*
        //testing
        var id1 = angular.element(document).find("ng-app-clothoWidgets");
        var id2 = id1.append(angular.element('<div clotho-widget clotho-widget-uuid="'+appUUID+'" clotho-widget-name="'+appInfo.moduleName+'"></div>'));
        var id3 = id2.append('<div ng-view></div>');

        console.log(id1);
        console.log(id2);
        console.log(id3);
        */

        //angular version
        //note angular returns parent, not appended element
        //todo - if want this, select appropriate child element
        //var insertInto = angular.element(document).find("ng-app-clothoWidgets").append(angular.element('<div clotho-widget clotho-widget-uuid="'+appUUID+'" clotho-widget-name="'+appInfo.moduleName+'"></div>').append('<div ng-view></div>'));

        //jQuery version
        var insertInto = $($('<div clotho-widget clotho-widget-uuid="'+appUUID+'" clotho-widget-name="'+appInfo.moduleName+'"></div>').append('<div ng-view></div>')).appendTo($clotho.appWidgets);

        $.getScript(appInfo.moduleUrl)
            .done(function () {
                angular.bootstrap(insertInto, [appInfo.moduleName]);
                deferred.resolve([appUUID, "[clotho-widget-uuid="+appUUID+"]"]);
                $rootScope.$safeApply();
            });

        return deferred.promise;
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
        fn.api.emit('log', msg);
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
        fn.api.emit('say', packaged);
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
        fn.api.emit('notify', data);
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
        fn.api.emit('alert', packaged);
    };

    //testing - search bar commands
    //todo - move outside of Clotho API
    var autocomplete = function(query) {
        var packaged = {
            "query" : query
        };

        //todo - use angular $cacheFactory to cache searches

        fn.searchbar.emit('autocomplete', packaged);
    };
    var autocompleteDetail = function(uuid) {
        //check the collector first
        if (Collector.hasItem('detail_'+uuid)) {
            return Collector.retrieveModel('detail_'+uuid);
        }

        var deferred = $q.defer();
        var packaged = {
            "uuid" : uuid
        };
        fn.searchbar.emit('autocompleteDetail', packaged);

        //testing
        PubSub.once('autocompleteDetail_function_id123', function(data) {
        // PubSub.once('autocompleteDetail_'+uuid, function(data) {
            $rootScope.$safeApply(deferred.resolve(data));
        }, '$clotho');

        return deferred.promise;
    };

    var submit = function(query) {
        var packaged = {
            "query" : query
        };
        fn.searchbar.emit('submit', packaged);
    };

    /**
     * @name Clotho.run
     *
     * @param {string} assistant    Function to be run on the server
     * @param {object} data JSON passed to the function
     * @param {boolean=} async   True if the return is asynchronous or false for synchronous (default)
     * @param {function()=} callback   Function to run upon promise resolving/failing
     *
     * @description
     * Runs a function on the server with the passed data
     *
     * @return {object} result of running the function on the server. Returns a promise if `async` is true, otherwise waits till completion to return.
     */
    var run = function clothoAPI_run(assistant, data, async, callback) {
        async = async || true;
        var packaged = {
            "assistant" : assistant || "",
            "data" : msg,
            "async" : async
        };
        fn.api.emit('run', packaged);
    };

    /**
     @name Clotho.recent
     *
     * @description
     * Request your most recently / commonly used sharables
     */
    var recent = function() {
        fn.api.emit('requestRecent', {});

        var deferred = $q.defer();

        PubSub.once('displayRecent', function(data) {
            $rootScope.$safeApply(deferred.resolve(data));
        }, 'clothoAPI');

        return deferred.promise;
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
        fn.api.emit('learn', packaged);
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
     * @param {object} quiz
     *
     * @description
     * Pass results of a quiz to the server to be graded
     *
     */
    var gradeQuiz = function clothoAPI_gradeQuiz(quiz) {
        fn.api.emit('gradeQuiz', quiz);

        var deferred = $q.defer();

        PubSub.once('quizResult:' + quiz.uuid, function(data) {
            $rootScope.$safeApply(deferred.resolve(data));
        }, 'clothoAPI');

        return deferred.promise;
    };


    return {
        //api
        get : get,
        get_url : get_url, //testing
        get_template : get_template, //testing
        get_script : get_script, //testing
        set : set,
        //clone : clone,
        query : query,
        create : create,
        edit : edit,
        revert : revert,
        destroy : destroy,
        show : show,
        show_simple : show_simple, //testing
        say : say,
        log : log,
        alert : alert,
        run : run,
        recent: recent,
        notify : notify,
        gradeQuiz : gradeQuiz,

        //toolkit
        bootstrap: bootstrap,
        watch : watch,
        watch2 : watch2,
        listen : listen,
        silence : silence,
        trigger: trigger,
        emit : emit,
        broadcast : broadcast,
        on : on,
        once : once,
        off : off,

        //searchbar
        submit: submit,
        autocomplete : autocomplete,
        autocompleteDetail : autocompleteDetail

    }
}]);

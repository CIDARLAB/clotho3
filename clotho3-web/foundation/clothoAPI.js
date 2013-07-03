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
 *  Most functions are written so that resulting behavior is not captured in a return, but rather by behavior within the ClientAPI. For some cases, returns are necessary for the function to make sense (e.g. Clotho.get() ).
 *
 *  At present, the goal is write them such that all resulting behavior will be captured by the messages returned over the socket from the server -- i.e. behaviour within clientAPI, PubSub, etc. This may not be possible, but regardless, the channel '$clotho' is reserved for the serverAPI (and falls through on clientAPI).
 * The main problem is get()... which by its nature should return the object to the callback. However, it may make sense to encapsulate this functionality in the collector. See note below
 *
 * Future todos
 *  - return angular.noop for empty callbacks, or deal with better
 */

//todo - incorporate requestId field

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
     * @param {string|array} uuid UUID of Sharable, or an array of string UUIDs
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
     * @param {object} sharable  JSON representation of Sharable containing new data - may be partial (just a few fields), but must contain an ID or UUID
     *
     * @description
     * Updates a sharable with passed object. Sets the fields present in the spec to the values in the spec. Fields missing from the spec are unchanged. If a field is missing from the specified object, it will be created. Clotho will send an error message on the 'say' channel if an update could not be applied.
     *
     * Server will emit a collect(uuid) call, upon object being updated
     */
    var set = function clothoAPI_set(sharable) {
        //strip $$v in case use promise to set multiple times
        while (sharable.$$v) {
            sharable = sharable.$$v;
        }

        fn.api.emit('set', sharable);
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



    /**
     * @name Clotho.query
     *
     * @param obj
     *
     * @description
     * Returns all objects that match the fields provided in the spec.
     */
    var query = function(obj) {
        fn.searchbar.emit("query", obj);
        var deferred = $q.defer();

        PubSub.once('query', function(data) {
            $rootScope.$safeApply(deferred.resolve(data));
        }, '$clotho');

        return deferred.promise;
    };
    /**
     * @name Clotho.create
     *
     * @param {object} object A json object or list of json objects that describe(s) an instance or instances of a clotho schema.
     *
     * @description
     * Creates the described objects. Returns a list of uuids of the created objects, in the same order as they were sent in. If an object was not created, that object's entry in the list of uuids will be null. Reasons for failing to create the failed objects will be provided in the 'say' channel.
     *
     * @returns {object} The object if created successfully, null otherwise
     */
    var create = function clothoAPI_create(object) {

        fn.api.emit('create', object);

        var deferred = $q.defer();

        PubSub.once('create', function(data) {
            $rootScope.$safeApply(deferred.resolve(data));
        }, '$clotho');

        return deferred.promise;
    };

    /**
     * @name Clotho.destroy
     *
     * @param {string|array} uuid UUID of Sharable, or an array of string UUIDs
     *
     * @description
     * Destroys the listed objects. Does not destroy anything if the selector is ambiguous. Clotho will send an error message on the 'say' channel if destroy fails for an object.
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
     * @param {object} parameters Example given below
     *
     format:
     {
         "template" : <url>,         // required
         "target" : <DOM ELEMENT>    // suggested, or absolute positioning in CSS
         "args" : {<object>}         // data to copy onto $scope
         "controller" : <url>,       // optional
         "dependencies" : [
             <urls>                  // required if in controller
         ],
         styles : {
             <styles>
             //e.g.
             "background-color" : "#FF0000"
         }
     }

     note CAVEATS:
     - currently, controllers etc. must be tied to Application.Extensions.___

     */
    var show = function(parameters) {
        fn.api.emit('show', parameters);
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
        PubSub.once('autocompleteDetail_'+uuid, function(data) {
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
     * @param {string} func Object UUID or name indicating the function to be run (follows get semantics for ambiguous selectors)
     * @param {object} args A JSON object with key-value pairs providing the argument values to the function.
     *
     * @description
     * Runs a function on the server on supplied arguments.
     * Clotho will send an error message on the 'say' channel if there is an error during function execution, there is no function matching the specifier, or there are ambiguously specified arguments.
     */
    var run = function clothoAPI_run(func, args) {
        var packaged = {
            "function" : func || "",
            "arguments" : args
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
            console.log('runAPI');
            console.log(data);
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
        show_old : show_old,
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

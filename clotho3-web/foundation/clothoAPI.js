'use strict';

/**
 * @name serverAPI
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

Application.Foundation.service('Clotho', ['Socket', 'Collector', 'PubSub', '$q', '$rootScope', '$location', function(Socket, Collector, PubSub, $q, $rootScope, $location) {

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
     * @returns {Promise | object} Promise to be fulfilled with JSON representation of Sharable, or the object for synchronous (blocking)
     *
     * todo - add callback? Promises too complicated?
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
            fn.emit('get', uuid);

            PubSub.once('model_change:'+uuid, function(data) {
                //console.log("CLOTHOAPI\t");
                //console.log(data);
                $rootScope.$apply(deferred.resolve(data));
            }, '$clotho');

            if (synchronous) {
                //fixme need to timeout or something
                return deferred.promise.then(function(result) {
                    return result;
                });
            } else {
                return deferred.promise;
            }
        }
    };

    //for testing
    var get_url = function clothoAPI_get_url(uuid) {
        var deferred = $q.defer();

        fn.emit('get:url', uuid);

        PubSub.once('model_change:'+uuid, function(data) {
            $rootScope.$apply(deferred.resolve(data));
        });

        return deferred.promise;
    };

    //for testing
    var get_template = function clothoAPI_get_template(uuid) {
        fn.emit('get:template', uuid);
    };

    //for testing
    var get_script = function clothoAPI_get_script(uuid) {
        fn.emit('get:script', uuid);
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
        //todo - do this natively, this is a hack
        //strip $$v in case use promise to set multiple times
        while (data.$$v) {
            data = data.$$v;
        }

        var packaged = {
            "uuid" : uuid,
            "data" : data
        };
        fn.emit('set', packaged);
    };

    var query = function(str) {

    };

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
        }, reference)
    };

    /**
     * @name Clotho.listen
     *
     * @param {string} channel PubSub Channel to listen to
     * @param {function} callback Function to be executed on change
     * @param {string=} reference Reference (e.g. $scope.$id) for element to unlink listener on destroy
     *
     * @description
     * Watches for published events on a given channel, and runs a callback on the event
     *
     * todo - smarter reference logic
     */
    var listen = function clothoAPI_listen(channel, callback, reference) {
        reference = typeof reference != 'undefined' ? reference : null;
        PubSub.on(channel, function clothoAPI_listen_callback(data) {
            $rootScope.$safeApply(callback(data));
        }, reference)
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
        fn.emit(channel, args);
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
        fn.emit('create', packaged);
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
        fn.emit('destroy', uuid);
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
        //fn.emit('edit', uuid);

        $location.path("/editor/" + uuid);
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
        fn.emit('show', packaged);
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
     */
    var nextWidgetUUID = 0;
    var bootstrap = function clothoAPI_bootstrap(appInfo) {

        nextWidgetUUID++;
        var appUUID = nextWidgetUUID;

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
        var insertInto = $($('<div clotho-widget clotho-widget-uuid="'+appUUID+'" clotho-widget-name="'+appInfo.moduleName+'"></div>').append('<div ng-view></div>')).appendTo('ng-app-clothoWidgets');

        $script(appInfo.moduleUrl, function() {
            angular.bootstrap(insertInto, [appInfo.moduleName]);
        });

        return [appUUID, "[clotho-widget-uuid="+appUUID+"]"];
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

    //todo - move outside of Clotho API
    var autocomplete = function(query) {
        var packaged = {
            "query" : query
        };
        fn.emit('autocomplete', packaged);
    };
    var autocompleteDetail = function(uuid) {
        var packaged = {
            "uuid" : uuid
        };
        fn.emit('autocompleteDetail', packaged);
    };
    var submit = function(query) {
        var packaged = {
            "query" : query
        };
        fn.emit('submit', packaged);
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
            "assistant" : userID,
            "data" : msg,
            "async" : async
        };
        fn.emit('run', packaged);
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
        fn.emit('learn', packaged);
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

    };

    /**
     * @name Clotho.gradeQuiz
     *
     * @param {string} uuid
     * @param {object} obj
     * @param {function} callback
     *
     * @description
     * Pass results of a quiz to the server to be graded
     *
     */
    var gradeQuiz = function clothoAPI_gradeQuiz(uuid, obj, callback) {

    };


    return {
        get : get,
        get_url : get_url, //testing
        get_template : get_template, //testing
        get_script : get_script, //testing
        set : set,
        query : query,
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
        create : create,
        edit : edit,
        destroy : destroy,
        show : show,
        bootstrap: bootstrap,
        say : say,
        log : log,
        alert : alert,
        autocomplete : autocomplete,
        autocompleteDetail : autocompleteDetail,
        submit: submit,
        run : run
    }
}]);
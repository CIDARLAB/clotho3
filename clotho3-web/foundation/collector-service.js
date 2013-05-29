'use strict';

//note : localStorage only supports strings, so we need to manually serialize and deserialize


Application.Foundation.service('Collector', ['$window', 'PubSub',function($window, PubSub) {

    return (window.$clotho.$collector) ? window.$clotho.$collector : window.$clotho.$collector = generateCollector();

    function generateCollector() {

        /************
         LOCAL STORAGE INTERFACE
        ************/

        var angularLS = generateAngularLS();

        function generateAngularLS() {

            //defaults
            var refStorage = $window.localStorage;
            var serializer = JSON;
            //prefix to prevent name-clashes
            var prefix = "clotho_";


            // --- Check for support ---
            // FUTURE - add this in later to be more robust
            // e.g. https://github.com/grevory/angular-local-storage/blob/master/localStorageModule.js
            // and could fallback to cookies for small strings (<4kb)

            // Checks the browser to see if local storage is supported
            var browserSupportsLocalStorage = function () {
                try {
                    return ('localStorage' in window && window['localStorage'] !== null);
                } catch (e) {
                    return false;
                }
            };

            // Checks the browser to see if cookies are supported
            // look into angular $cookies if we want to implement this as a backup
            // note the limitations of cookies before committing to that
            var browserSupportsCookies = function() {
                try {
                    return navigator.cookieEnabled ||
                        ("cookie" in document && (document.cookie.length > 0 ||
                            (document.cookie = "test").indexOf.call(document.cookie, "test") > -1));
                } catch (e) {
                    return false;
                }
            };

            //testing
            //console.log("localStorage support? " + browserSupportsLocalStorage());

            // --- local storage interface ----

            var clear = function() {
                refStorage.clear();
                return( this );
            };

            // returns an item, or optional defaultValue if not found
            var getItem = function( key, defaultValue ){
                var value = refStorage.getItem( prefix+key );
                if (value == null){
                    // check to see if default is falsy
                    return((typeof( defaultValue ) != "undefined") ? defaultValue : null );
                } else {
                    return(serializer.parse( value ) );

                }
            };

            // return a boolean whether a given object exists
            var hasItem = function( key ){
                return( refStorage.getItem( prefix+key ) != null );
            };

            //remove a given item from storage
            var removeItem = function( key ){
                refStorage.removeItem( prefix+key );
                return( this );
            };

            //adds an item to storage, automatically serializing.
            //NOTE -  cannot add functions or private variables
            var setItem = function( key, value ){
                //prevent adding of empty objects
                if (!value && value!==0 && value!=="") return false;
                refStorage.setItem(prefix+key, serializer.stringify( value ) );
                return( this );
            };

            // ----- LOCAL STORAGE EVENTS ----
            /* NOTES
             this half is for model update broadcasting across pages (via localStorage)
             the other half is using PUB/SUB (within pages + to bubble up updates)

             note - implementation of localStorage events varies & is unreliable across browsers
             - event handlers only invoked for current window / tab where data is written / deleted (even though spec states all windows)
             - may not be called when a key is updated, and only when it is created / deleted
             - storage deletion only returns key, not deleted value
             */

            var handle_storage_change = function(e) {
                //console.log("change made to local storage");
                if (!e) { e = window.event; }

                //TODO - better checking for e.key
                var uuid = e.key.replace(prefix, '') || '';
                console.log("handle_storage_event");
                PubSub.trigger("model_change", uuid);
                PubSub.trigger("model_change:"+uuid, getItem(uuid));
            };

            if (window.addEventListener) {
                window.addEventListener("storage", handle_storage_change, false);
                //testing
                // console.log("local storage Listener added");
            } else {
                window.attachEvent("onstorage", handle_storage_change); //IE8
                //testing
                // console.log("local storage Listener added for IE8");
            }


            return {
                getPrefix : function() {return prefix},
                isSupported : browserSupportsLocalStorage,
                clear : clear,
                getItem : getItem,
                hasItem : hasItem,
                removeItem : removeItem,
                setItem : setItem
            }
        }

        /*******
        COLLECTOR
        *******/

        // ------- DATA STORAGE + ACCESS -----

        // collector is our internal model / hashmap
        // localStorage only supports strings, so must be stored separately (can't use same obj reference)
        var collector = {};

        //does not broadcast update (locally at least, but angularLS will via localStorage updates)
        var silentAddModel = function(uuid, obj) {
            collector[uuid] = obj;
            angularLS.setItem(uuid, obj);
        };

        //passes update message - usual way of adding model to collector
        var storeModel = function(uuid, obj) {
            //testing
            //console.log("COLLECTOR\tstoring uuid " + uuid);
            //console.log(obj);
            silentAddModel(uuid, obj);
            PubSub.trigger("model_change", uuid);
            PubSub.trigger("model_change:" + uuid, obj);
        };

        //check if we have a model. Return it if present. Otherwise, return false. Returns the collector if no UUID is passed.
        var retrieveModel = function(uuid) {
            if (typeof uuid == 'undefined') {
                return collector;
            }
            if (collector[uuid] && typeof collector[uuid] != 'undefined')
                return collector[uuid];

            if ( collector[uuid] = angularLS.getItem(uuid) ) {
                return collector[uuid];

            } else {
                return false;
            }
        };

        var clearStorage = function() {
            angularLS.clear();
            collector = {};
            PubSub.trigger("collector_reset", null);
        };

        /*
         future - may want to implement something like this once websocket is set up
         e.g. sync with server... this is constant, may need to throttle if too much communication

         fail-safe in case PubSub doens't pick something up
         note - watching collector, so localStorage changes not picked up here (will be in angularLS)

         $rootScope.$watch('collector', function() {

         });
         */

        // ---- PUB / SUB ----

        //todo - need to rewrite (not correct right now), move to clotho API
        /**
         * @name Collector.watch
         *
         * @param toWatch {string | object} event to watch.
         *  - If pass `all`, then all model changes will be passed to toUpdate
         *  - Otherwise, defaults to the UUID you pass in
         * @param callback {function} function to be run. Passed new model.
         *
         * @description
         * Function to be inherited by controllers / services, which will run a given callback upon the changing of the given uuid.
         */
        var watch = function(toWatch, callback) {
            //handle watching all model changes
            if (toWatch == "all") {toWatch = "model_change"}
            else {toWatch = "model_change:" + toWatch}

            PubSub.on(toWatch, function(uuid) {
                console.log("COLLECTOR.update\tuuid: " + toWatch);
                callback(retrieveModel(uuid));
            });

        };

        // ------- FACADE -----------

        return {
            collector : collector,
            silentAddModel : silentAddModel,
            storeModel : storeModel,
            retrieveModel : retrieveModel,
            clearStorage : clearStorage,
            watch : watch
        }
    }
}]);
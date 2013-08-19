'use strict';

Application.Services.service('CollectorOld', ['$resource', '$rootScope', 'angularLS', 'PubSub', '$q', function($resource, $rootScope, angularLS, PubSub, $q) {
    // ----- TODO ----
    // way to request arrays (instead of only objects)
    // avoid UI flickers
    // FUTURE
    // check out cacheFactory later - angular support for something similar? docs are minimal right now, maybe built out more later...

    // note - local storage done behind the scenes, Collector checks there first, otherwise requests remote file

    console.log("localStorage support? " + angularLS.isSupported());

    // ------- DATA STORAGE + ACCESS -----

    // collector is our internal model / hashmap
    // localStorage only supports strings, so must be stored separately (can't use same obj reference)

    var collector = {};

    //does not broadcast update (locally at least, but angularLS will via localStorage updates)
    var silentAddModel = function(uuid, obj) {
        if (!uuid) return;
        if (typeof obj == 'string') {obj = JSON.parse(obj);}
        collector[uuid] = obj;
        angularLS.setItem(uuid, obj);
    };

    //passes update message - usual way of adding model to collector
    var storeModel = function(uuid, obj) {
        console.log("COLLECTOR\tstoring uuid " + uuid);
        console.log(obj);
        silentAddModel(uuid, obj);
        PubSub.publish("model_change", uuid);
        PubSub.publish("model_change:" + uuid, obj);
    };

    //temporary - necessary right now for retrieving model (not yet a separate service)
    //when we move to more robust resource requests, look into $cacheFactory
    // $cacheFactory is a built in Angular cache used by $http... minimal docs right now though
    var RemoteModel = $resource('models/:unique_id.json', {}, {
        //getFirst: {method:'GET', params:{uuid:'first'}, isArray:false}
    });

    //for requesting any URL
    //e.g. SimpleReq.get({url: <url>})
    //future - implement safely -- probably need to do when have real server
    var SimpleReq = $resource(':url', {}, {});

    //returns a promise
    var remoteModel = function(uuid) {
        var deferred = $q.defer();

        RemoteModel.get({unique_id:uuid}, function(data) {
            //if you want to modify the data or something... do it here
            storeModel(uuid, data);
            deferred.resolve(data);
            //updating the view is handled in storeModel (broadcast published) - currently requires subscription to model update

        }, function() {
            deferred.reject();
        });
        return deferred.promise;
    };

    // get a model, checking local storage first, otherwise requesting it
    var retrieveModel = function(uuid) {
        var deferred = $q.defer();
        var promise = deferred.promise;

        if (typeof uuid === 'undefined') {
            console.log("retrieveModel called without param");
            deferred.resolve( collector );
        }

        if ( collector[uuid] = angularLS.getItem(uuid) ) {
            //object was in local storage, return the data
            promise.then(function() {
                //testing - fn() will run... if wrap in function, (i.e. function(){.code.} ) doesn't run until resolved
                //deferred.resolve( collector[uuid] )
                console.log("COLLECTOR\tpromise test\tthis is run --- " + uuid)
            });
            deferred.resolve( collector[uuid] );
            //testing:
            console.log("COLLECTOR\tmodel in local storage: " + uuid);
            console.log(deferred);

        } else {
            //doesn't exist, request it
            remoteModel(uuid).then(function (val) {
                deferred.resolve(val);
            });
            //testing:
            console.log("COLLECTOR\tmodel is remote model: " + uuid);
            console.log(deferred);
        }
        //testing
        console.log("COLLECTOR\tretrieve model promise: " + uuid);
        console.log(deferred);

        //note - as currently exists, need to use a then() statement on the promise, or may not work properly. Don't really understand. figure this out.
        //fixme -- want to return promise, but if resolved before return (e.g. in localStorage) then will pass the promise, and it won't be resolved after the return
        // see SO article: http://stackoverflow.com/questions/15802553
        return promise;
    };

    //future - implement. should have base path so can't download arbitrary files
    //download a passed URL and store for a given UUID
    var downloadModel = function(uuid, url) {
        console.log("COLLECTOR DOWNLOAD url: " + url);
        var deferred = $q.defer();
        SimpleReq.get({url: url}, function(data) {
            storeModel(uuid, data);
            deferred.resolve(data);
        }, function() {
            deferred.reject();
        });
        return deferred.promise;
    };


    var clearStorage = function() {
        angularLS.clear();
        PubSub.publish("collector_reset", null);
    };

    // future - implement something like this once websocket is set up
    $rootScope.$watch('collector', function() {
        //sync with server... this is constant, may need to throttle if too much communication

    });

    // ---- PUB / SUB ----

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

        PubSub.subscribe(toWatch, function(uuid) {
            console.log("COLLECTOR.update\tuuid: " + toWatch);
            retrieveModel(uuid).then(function(result) {
                callback(result);
            });
        });

    };

    // ------- FACADE -----------

    return {
        collector : collector,
        silentAddModel : silentAddModel,
        storeModel : storeModel,
        resetModel : remoteModel,
        retrieveModel : retrieveModel,
        //downloadModel : downloadModel,
        clearStorage : clearStorage,
        watch : watch
    }

}]);
angular.module('clotho.core').service('Collector',
	function($window, clothoLocalStorage, PubSub, Debug) {

    return ($window.$clotho.$collector) ? $window.$clotho.$collector : $window.$clotho.$collector = generateCollector();

    function generateCollector() {

        function broadcastModelUpdate(uuid, obj) {
            PubSub.trigger("update", uuid);
            PubSub.trigger("update:" + uuid, angular.copy(obj));
        }

	      var Debugger = new Debug('Collector', '#55bb55');

        /*******
        COLLECTOR
        *******/

        // ------- DATA STORAGE + ACCESS -----

        //does not broadcast update (locally at least, but angularLS will via localStorage updates)
        var silentAddModel = function(uuid, obj) {
            clothoLocalStorage.setItem(uuid, obj);
        };

        //passes update message - usual way of adding model to collector
        //pass true for 'force' to force collect even if obj identical and broadcast of update
        var storeModel = function(uuid, obj, force) {
	          //todo - ensure that what is in collector also matches localStorage
            if (force || !angular.equals(retrieveRef(uuid), obj)) {
	            Debugger.log(uuid + " (saving)", obj);
                silentAddModel(uuid, obj);
                broadcastModelUpdate(uuid, obj);
            }
            else {
	            Debugger.log(uuid + " (model unchanged)");
            }
        };

        //check if we have a model. Return it if present. Otherwise, return false. Returns the collector if no UUID is passed.
        var retrieveModel = function(uuid) {
            return angular.copy(retrieveRef(uuid));
        };

        var retrieveRef = function(uuid) {
            if ( clothoLocalStorage.hasItem(uuid)) {
	            return clothoLocalStorage.getItem(uuid)
            } else {
                return false;
            }
        };

        var removeModel = function (uuid) {
            if (clothoLocalStorage.hasItem(uuid)) {
                clothoLocalStorage.removeItem(uuid);
                return true;
            } else {
                return false;
            }
        };

        var clearStorage = function() {
            clothoLocalStorage.clear();
            PubSub.trigger("collector_reset");
        };

	      //todo - need a way to only clear if this is the first tab of clotho opening -- see #310
	      //note - clear on first tab because want a fresh instance of clotho - don't want to load stale models
	      clearStorage();


        // ------- FACADE -----------

        return {
            hasItem : clothoLocalStorage.hasItem,
            silentAddModel : silentAddModel,
            storeModel : storeModel,
            retrieveModel : retrieveModel,
            retrieveRef : retrieveRef,
            removeModel : removeModel,
            clearStorage : clearStorage
        }
    }
});
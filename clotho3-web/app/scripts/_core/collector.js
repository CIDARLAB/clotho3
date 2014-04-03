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

        // collector is our internal model / hashmap
        // localStorage only supports strings, so must be stored separately (can't use same obj reference)
        var collector = {};

        //does not broadcast update (locally at least, but angularLS will via localStorage updates)
        var silentAddModel = function(uuid, obj) {
            collector[uuid] = obj;
            clothoLocalStorage.setItem(uuid, obj);
        };

        //passes update message - usual way of adding model to collector
        //pass true for 'force' to force collect even if obj identical and broadcast of update
        var storeModel = function(uuid, obj, force) {
	          //todo - ensure that what is in collector also matches localStorage
            if (force || !angular.equals(retrieveRef(uuid), obj)) {
	            Debugger.log(uuid + " (saving)", collector[uuid], obj);
                silentAddModel(uuid, obj);
                broadcastModelUpdate(uuid, obj);
                //testing Debugger.log(collector[uuid]);
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
            if (typeof uuid == 'undefined') {
                return collector;
            }
            if (collector[uuid] && typeof collector[uuid] != 'undefined')
                return collector[uuid];

            if ( collector[uuid] = clothoLocalStorage.getItem(uuid) ) {
                return collector[uuid];

            } else {
                return false;
            }
        };

        var removeModel = function (uuid) {
            if (collector[uuid]) {
                collector[uuid] = null;
                clothoLocalStorage.removeItem(uuid);
                return true;
            } else {
                return false;
            }
        };

        var clearStorage = function() {
            clothoLocalStorage.clear();
            collector = {};
            PubSub.trigger("collector_reset");
        };


        // ------- FACADE -----------

        return {
            collector : collector,
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
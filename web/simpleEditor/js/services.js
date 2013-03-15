'use strict';

/* Services */

angular.module('editorServices', ['ngResource'])
    .factory('Institution', function($resource, $rootScope) {

        // ------- DATA ---------

        //just an example JSON... realistically we'll probably always be pulling from the server.
        //as an aside, could merge all these address fields into one if matches server model better
        var localExample = {
            "displayName" : "John Smith",
            "city" : "Berkeley",
            "state" : "CA",
            "country" : "USA",
            "noAccessParam" : "Can't touch this"
        };

        //the model for the views, shared across controllers
        //could store as a hash, but once we have a server that allows us to actually save data it'd be pointless
        var demoData;

        //which institutionID the user is using
        var demoInstID;

        // ------- AJAX ---------

        //could be a separate service, but I want to encapsulate the logic in here for now - no need to separate
        //
        //getFirst() - custom method with hardcoded instID as example
        //
        //we can define more methods like this - using clotho jargon ...... BUT
        //the main issue is it's difficult to know which controller OR service to update,
        //and it requires much more custom code / routing to go in/out of Angular so much
        //
        //Right now, I think it's easier to stick to $resource... header functionality should be more apparent in 1.1.0
        //plus its nice that it handles creating a factory for us. Can change to $http() if necessary
        var RemoteInstitution = $resource('institutions/:instID.json', {}, {
            getFirst: {method:'GET', params:{instID:'first'}, isArray:false}
        });

        //return the model data (not a reference to it) - not sure if this is the best architecture though...
        var getModel = function(institutionID) {
            //demonstrate $scope lifecycle - this will get called twice if you uncomment
            //console.log("demoInstID : " + demoInstID + "\tinstitutionID : " + institutionID);

            //first, check if we're requesting a new institution ID (or we haven't requested one yet)
            if (demoInstID != institutionID || typeof demoInstID === 'undefined') {
                demoInstID = (typeof institutionID === 'undefined') ? 'first' : institutionID;
                getRemote(demoInstID);
            }
            //if not requesting something new, then pass demoData via broadcast
            //later, not really an issue because want to pull info from server,
            //but right now want to persist across controllers
            else {
                sendBroadcast(MODEL_CHANGE);
            }

        };

        //deals with promise callback
        //institutionID e.g. using $routeParams to get a specific one
        var getRemote = function(institutionID) {
            if (typeof institutionID === 'undefined') {
                console.log("getRemote called without param");
                return;
            }
            RemoteInstitution.get({instID:demoInstID}, function(data) {
                //if you want to modify the data or something... do it here
                demoData = data;
                sendBroadcast(MODEL_CHANGE);
            });
        };

        //changing institutions without location change
        var setInst = function(institutionID) {
            if (typeof institutionID !== 'undefined' && demoInstID != institutionID) {
                demoInstID = institutionID;
                getRemote(demoInstID);
            }
        };

        // -------- PUB / SUB -------

        //local constants, private
        var MODEL_CHANGE = 'modelChange';

        //PUBLISH
        //sends a passed message to controllers, which do something, like update their model...
        //
        //exposed to controllers indirectly, e.g. in save()
        //
        //this is sort of a master controller, which can broadcast multiple messages
        //
        //Angular requires an explicit update (either like this, or using $apply() or $digest())
        //when something is done outside of Angular....
        //main examples are AJAX, DOM Manipulation, timeouts
        var sendBroadcast = function(subject) {
            $rootScope.$broadcast(subject);
        };

        //not currently in use
        //pass in an object to set as demoData, broadcast change
        //may be useful if we need more logic than just persisted in the model (currently not using)
        var prepBroadcast = function(channel, instObj) {
            if (typeof instObj !== 'undefined') {demoData = instObj;}
            sendBroadcast(channel);
        };

        //SUBSCRIBE
        //pass a subscribe method to controllers (included in facade)
        //purpose: check for state persistence across controllers, and...
        //require $scope to avoid creation of closure when destroyed
        var onDataChange = function($scope, handler) {
            $scope.$on(MODEL_CHANGE, function(){
                // pass the data to the handler
                handler(angular.copy(demoData));
            });
        };

        // ----- FACADE -------

        //methods exposed to controllers
        return {
            getModel : getModel,

            getRemote : getRemote,

            setInst : setInst,

            reset : function() {
                //as example, could be remote which probably makes more sense once actually serving data
                //i.e. so make call to server for object each time want to reset
                demoData = getRemote(demoInstID);
                //internal publish
                sendBroadcast(MODEL_CHANGE);
            },

            discard : function(newData) {
                console.log("discarding! discarded object follows: ");
                console.log(newData);
            },

            save : function(newData) {
                demoData = newData;

                console.log("saving! Object follows: ");
                console.log(demoData);

                // this is where we would send something to the server
                // ...i.e. where would probably POST or PUT something



                //internal publish
                sendBroadcast(MODEL_CHANGE)
            },

            //expose subscribe
            onDataChange : onDataChange
        }
    });
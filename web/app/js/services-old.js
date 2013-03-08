'use strict';

/* Services */

angular.module('editorServicesOLD', ['ngResource'])
    .factory('Institution', function($resource) {
        //could merge all these address fields into one if matches server model better
        var localMaster = {
            "displayName" : "John Smith",
            "city" : "Berkeley",
            "state" : "CA",
            "country" : "USA",
            "noAccessParam" : "Can't touch this"
        };

        var RemoteInstitution = $resource('institutions/:instID.json', {}, {
            query: {method:'GET', params:{instID:'first'}, isArray:false}
        });

        var demoData;

        //example of asynchronous call - populate data with callback on request return
        RemoteInstitution.query({}, function(data) {
            demoData = data;
            return demoData;
        });

        return {
            getInst : function() {
                return demoData;
            },

            getLocal : function() {
                return angular.copy(localMaster);
            },

            getRemote : function() {
                //change e.g. using $routeParams to get a specific one
                return RemoteInstitution.query();
            },

            reset : function() {
                //as example, could be remote which probably makes more sense once actually serving data
                demoData = angular.copy(localMaster);
                return demoData;
            },

            save : function() {
                var temp = angular.copy(demoData);
                //do something meaningful here
                console.log(temp);
            }
        }
    });

'use strict';

Application.Extensions.controller('writingTrailsCtrl', ['$compile', '$scope', 'Clotho', '$http', function($compile, $scope, Clotho, $http) {
    $scope.getSchema = function(schemaName) {
        return Clotho.query({schema: "ClothoSchema", name: schemaName})
            .then(function(results) {
                console.log(results[0]);
                return results[0]
            });
    };

    $scope.getTrail = function(trailName) {
        return Clotho.query({schema: 'Trail', name: trailName})
            .then(function(results) {
                console.log(results[0]);
                return results[0]
            });
    };

    //temporary hack while Trail schema can't be imported
    $scope.getTrailHack = function() {
        return $http.get('/models/definition_trail.json').then(function(data) {
            return data.data;
        });
    }


}]);
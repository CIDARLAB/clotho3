'use strict';

Application.Extensions.controller('writingTrailsCtrl', ['$compile', '$scope', 'Clotho', function($compile, $scope, Clotho) {
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



}]);
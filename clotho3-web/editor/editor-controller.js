'use strict';

Application.Editor.controller('EditorCtrl', ['$scope', '$routeParams', '$location', 'Clotho', function($scope, $routeParams, $location, Clotho) {

    $scope.setUUID = function (uuid) {
        $scope.uuid = uuid;
        $scope.sharable = Clotho.get($scope.uuid);
        $location.path('/editor/' + $scope.uuid);
    };

    //init()
    $scope.setUUID($routeParams.uuid);


    //testing
    $scope.setNewFirst = function() {
        var tempObj = {
            "$clotho" : {
                "schema" : "schema_institution",
                "uuid" : "inst_first"
            },
            "displayName" : "New Person"
        };
        Clotho.set('inst_first', tempObj);
    };




}]);
'use strict';

Application.Editor.controller('EditorCtrl', ['$scope', '$routeParams', '$location', 'Clotho', function($scope, $routeParams, $location, Clotho) {

    $scope.setID = function (id) {
        $scope.id = id;
        $scope.sharable = Clotho.get($scope.id);
        $location.path('/editor/' + $scope.id);
    };

    //init()
    $scope.setID($routeParams.id);

}]);
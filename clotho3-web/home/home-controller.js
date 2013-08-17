'use strict';

Application.Primary.controller('HomeCtrl', ['$scope', '$location', function($scope, $location) {

    $scope.enterClotho = function() {
        $location.path('/trails/bb7f191e810c19729de860bb');
    };

    $scope.enterEugene = function() {
        $location.path('/trails/bb02191e810c19729de860aa');
    };

    $scope.enterRaven = function() {
        $location.path('/trails/bb02191e810c19729de860bb');
    };

}]);
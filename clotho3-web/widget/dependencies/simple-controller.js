'use strict';

Application.Widgets.controller('SimpleCtrl', ['$scope', 'SimpleService', function($scope, SimpleService) {
    $scope.test = "string works!";
    $scope.serviceText = SimpleService.text();
}]);
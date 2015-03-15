'use strict';

Application.Extensions.controller('SimpleCtrl', ['$scope', 'SimpleService', function($scope, SimpleService) {
    $scope.test = "string works!";
    $scope.serviceText = SimpleService.text();
}]);
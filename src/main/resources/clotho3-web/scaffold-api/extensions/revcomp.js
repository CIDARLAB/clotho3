'use strict';

Application.Extensions.controller('revcompCtrl', ['$scope', 'DNA', function($scope, DNA) {
    $scope.sequence = 'aaaaaaaggggggcccctttttt';
    $scope.DNA = DNA;
}]);
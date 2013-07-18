'use strict';

Application.Extensions.controller('revcompCtrl', ['$scope', 'DNA', function($scope, DNA) {
    $scope.DNA = DNA;
    $scope.sequence = 'acgtacgatcgatcgat';
}]);
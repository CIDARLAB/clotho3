'use strict';

console.log('revcomp Ctrl downloaded');

Application.Extensions.controller('revcompCtrl', ['$scope', 'DNA', function($scope, DNA) {

    console.log('revcompCtrl instantiated');

    $scope.DNA = DNA;
    $scope.sequence = 'acgtacgatcgatcgat';
}]);
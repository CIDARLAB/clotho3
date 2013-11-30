'use strict';

//todo - remove template reliance via $routeProvider -- just don't use ng-view
angular.module('revcompApp', ['clothoPackage'])
    .config(['$routeProvider', function ($routeProvider) {
        $routeProvider.
            otherwise({
                templateUrl:'partials/trails/trail_super/widget/default-partial.html'
            })
    }])
    .run(['$rootScope', 'Clotho', function($rootScope, Clotho) {
        $rootScope.Clotho = Clotho;
    }])
    .controller('revcompAppCtrl', ['$scope', 'Clotho', function($scope, Clotho) {
        $scope.sequence = "ttttttttttttttggggcccaaa";
    }])
    .filter('dnaUppercaseA', function() {
        return function(input) {
            return angular.lowercase(input).replace(/a/g, 'A');
        }
    });
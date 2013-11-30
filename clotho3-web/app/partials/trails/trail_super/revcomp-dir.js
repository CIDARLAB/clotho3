'use strict';

Application.Extensions.directive('dnaUppercaseRevcomp', ['DNA', function(DNA) {
    return {
        restrict: 'A',
        require: 'ngModel',
        link: function (scope, element, attrs, ngModelCtrl) {

            var fn = function(input) {
                return angular.uppercase(DNA.revcomp(input));
            };

            ngModelCtrl.$parsers.push(fn)
        }
    }
}]);
'use strict';

Application.Extensions.directive('simpleText', function() {
    return {
        restrict: 'EAC',
        link: function($scope, element, attrs, controller) {
            element.text('iambar');
        }
    };
});

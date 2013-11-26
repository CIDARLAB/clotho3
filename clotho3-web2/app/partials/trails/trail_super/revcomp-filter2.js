'use strict';

Application.Extensions.filter('dnaLowercaseG', function() {
    return function(input) {
        return angular.uppercase(input).replace(/G/g, 'g');
    }
});
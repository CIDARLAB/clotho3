'use strict';

Application.Extensions.filter('dnaUppercaseG', function() {
    return function(input) {
        return angular.lowercase(input).replace(/g/g, 'G');
    }
});
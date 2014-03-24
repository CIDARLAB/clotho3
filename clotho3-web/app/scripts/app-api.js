'use strict';

//Api-only stack
angular.module('clothoRoot', ['clotho.core']);
angular.element(document).ready(function() {
	angular.bootstrap(document, ['clothoRoot']);
});
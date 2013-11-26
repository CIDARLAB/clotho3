//angular components needed, add a couple methods to prototype
angular.module('clotho.angularAdditions', [
		'ngCookies', 'ngSanitize', 'ngRoute'
	])
	.config(function() {
		//angular function extensions
		var ext = {};
		ext.isEmpty = function(value) {
			return angular.isUndefined(value) || value === '' || value === null || value !== value;
		};
		ext.isScope = function(obj) {
			return obj && obj.$evalAsync && obj.$watch;
		};
		angular.extend(angular, ext);
	});
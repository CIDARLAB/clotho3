'use strict';

angular.module('clotho.webapp')
.controller('ExecutorCtrl', function ($scope, $route, $location, $filter, Clotho, ClothoSchemas) {

	/* route handling */

	//check for updates to URL - this will only be triggered when URL is not set by above action
	$scope.$on('$routeUpdate', function(scope, next, current) {
		var updateId = next.params.id;

		// only activate if new -- won't happen if dropdown changed this already
		if ($scope.functionId != updateId) {
			$scope.functionId = updateId;
		}
	});

	$scope.$on('$destroy', function () {
		$location.search('id', null).replace();
	});

	//initial query handling
	if ( $route.current.params.id ) {
		$scope.functionId = $route.current.params.id;
	}

	//check for updates to function ID, and retrieve it
	$scope.$watch('functionId', function (newval, oldval) {
		Clotho.get(newval).then(function(r) {
			if (ClothoSchemas.isFunction(r)) {
				$scope.function = r;
			}
		});
	});

	//functionality

	$scope.queryWrapper = function (text, schema) {
		return Clotho.autocomplete(text)
		.then(function (results) {

			if (angular.isUndefined(schema)) {
				return results;
			}

			return $filter('filter')(results, function (r) {
				if (schema == 'function') {
					return ClothoSchemas.isFunction(r);
				} else {
					return ClothoSchemas.isInstanceOfSchema(r, schema);
				}
			});
		});
	};

	$scope.onExecute = function (result) {

	};
});

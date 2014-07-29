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

	function resetFunctionArgs (functionObject) {
		functionObject = functionObject || $scope.function;
		$scope.functionArgs = {};
		/*angular.forEach(functionObject.args, function (arg) {
			$scope.functionArgs[arg.name] = null;
		});*/
	}

	//check for updates to function ID, and retrieve it
	$scope.$watch('functionId', function (newval, oldval) {
		Clotho.get(newval).then(function(r) {
			if (ClothoSchemas.isFunction(r)) {
				$scope.function = r;
				resetFunctionArgs(r);
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

	$scope.isPrimitiveField = ClothoSchemas.isPrimitiveField;
	$scope.schemaReadable = ClothoSchemas.mapSchemaIdToName;

	//future - filter
	$scope.capitalize = function (word) {
		return word.substring(0,1).toUpperCase() + word.substr(1);
	};

	function flattenArgs (functionObject, argsObject) {
		var arr = [];
		angular.forEach(functionObject.args, function (arg) {
			arr.push(argsObject[arg.name]);
		});
		return arr;
	}

	$scope.executeFunction = function () {
		Clotho.run($scope.functionId, flattenArgs($scope.function, $scope.functionArgs))
		.then(function (result) {
			$scope.functionResult = result;
		});
	};

	$scope.clearArguments = function () {
		resetFunctionArgs();
		$scope.functionResult = null;
	};

	$scope.setArgument = function (name, id) {
		Clotho.get(id, {mute : true}).then(function (r) {
			$scope.functionArgs[name] = r;
		});
	}

});

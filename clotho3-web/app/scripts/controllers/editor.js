'use strict';

angular.module('clotho.webapp').controller('EditorCtrl', function ($scope, $route, $location, Clotho, ClothoSchemas) {

	// todo - better handle dynamic routing
	// note - could use search param, and set reloadOnSearch to false
	/*$scope.$watch('editable', function (newval, oldval) {
	 $location.url('/editor/' + $scope.editable.id).replace();
	 });*/

	//init()
	//todo - test
	var queryResult = $route.current.locals.queryResult;
	queryResult && console.log('query result', queryResult);
	if (angular.isDefined(queryResult) && queryResult.length) {
		$scope.editable = queryResult[0];
	} else {
		$scope.editable = $route.current.params.id;
	}

	$scope.editModePass = false;

	//data

	$scope.objectTypes = ClothoSchemas.sharableTypes;

	$scope.schemas = [];
	ClothoSchemas.retrievedSchemas.then(function (schemas) {
		$scope.schemas = schemas;
	});

	//future - phase out (may be able to already)
	$scope.allInstances = [];
	Clotho.query({}).then(function (data) {
		_.each(data, function (item) {
			if (item.schema != "BuiltInSchema") {
				$scope.allInstances.push(item);
			}
		});
	});

	//functionality

	//todo - update pending #164 and #165, query over multiple fields, replace editableList
	$scope.queryNameWrapper = function(text) {
		return Clotho.query({name: text}).then(function (result) {
			return $filter('limitTo')(result, 10);
		})
	};

	$scope.createNew = function (type) {
		if (type == 'Instance') {
			$scope.chooseSubtype = !$scope.chooseSubtype;
		} else {
			$scope.editable = ClothoSchemas.createScaffold(type);
			$scope.chooseSubtype = false;
		}
	};

	$scope.createNewInstance = function (item, model, label) {
		$scope.editable = ClothoSchemas.createScaffold(model);
		$scope.chooseSubtype = false;
	};

	$scope.editExisting = function (item, model, label) {
		$scope.editable = item;
	};
});

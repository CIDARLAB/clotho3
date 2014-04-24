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
		console.log('\n\n\nEDITOR CONTROLLER', $scope.editable);
	} else {
		$scope.editable = $route.current.params.id;
		console.log('\n\n\nEDITOR CONTROLLER ID ', $scope.editable);
	}

	$scope.editModePass = false;

	//data

	$scope.objectTypes = ClothoSchemas.sharableTypes;

	$scope.schemas = [];
	ClothoSchemas.retrievedSchemas.then(function (schemas) {
		$scope.schemas = schemas;
	});

	//functionality

	//todo - move to autocomplete
	$scope.queryWrapper = function(text) {
		return Clotho.query({name: text}, {maxResults : 10}).then(function (results) {
			console.log(results);
			return results || [];
		});
	};

	$scope.createNew = function (type) {
		if (type == 'Instance') {
			$scope.chooseSubtype = !$scope.chooseSubtype;
		} else {
			$scope.editable = ClothoSchemas.createScaffold(type);
			$scope.chooseSubtype = false;
			$scope.editModePass = true;
		}
	};

	$scope.createNewInstance = function (item, model, label) {
		$scope.editable = ClothoSchemas.createScaffold(model);
		$scope.chooseSubtype = false;
		$scope.editModePass = true;
	};

	$scope.editExisting = function (item, model, label) {
		$scope.editable = item;
		$scope.editModePass = true;
	};
});

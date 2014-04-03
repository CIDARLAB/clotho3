'use strict';

angular.module('clotho.webapp').controller('EditorCtrl', function ($scope, $route, $location, Clotho) {

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

	$scope.objectTypes = {
		"Instance" : {
			readable : "Instance"
		},
		"Function": {
			readable : "Function",
			scaffold : {
				schema: "Function",
				language: "JSONSCHEMA"
			}
		},
		"Schema": {
			readable : "Schema",
			scaffold : {
				schema: "ClothoSchema",
				language: "JSONSCHEMA"
			}
		},
		"View" :  {
			readable : "View",
			scaffold : {
				schema: "View",
				language: "JSONSCHEMA"
			}
		}
	};

	//todo - move this to angular-global value, since used multiple places
	$scope.schemas = [];
	Clotho.query({"schema": "Schema"}).then(function (schemas) {
		$scope.schemas = schemas;
		_.remove($scope.schemas, function (schema) {
			return !!$scope.objectTypes[schema.name];
		});
	});

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

	function createForEditor (scaffold, type) {
		Clotho.create(scaffold)
		.then(function (id) {
			console.log('created object of type ' + type + ' with id: ' + id);
			$scope.editable = id;
			$scope.editModePass = true;
		});
	}

	$scope.createNew = function (type) {
		if (type == 'Instance') {
			$scope.chooseSubtype = !$scope.chooseSubtype;
		} else {
			createForEditor($scope.objectTypes[type].scaffold, type);
			$scope.chooseSubtype = false;
		}
	};

	$scope.createNewInstance = function (item, model, label) {
		createForEditor({schema: model, language: "JSONSCHEMA"}, label);
		$scope.chooseSubtype = false;
	};

	$scope.editExisting = function (item, model, label) {
		$scope.editable = item;
	};
});

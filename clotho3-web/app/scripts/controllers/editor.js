'use strict';

angular.module('clotho.webapp').controller('EditorCtrl', function ($scope, $routeParams, $location, Clotho) {

	//note - could use search param, and set reloadOnSearch to false
	/*$scope.$watch('editable', function (newval, oldval) {
	 $location.url('/editor/' + $scope.editable.id).replace();
	 });*/

	//init()
	$scope.editable = $routeParams.id;
	$scope.editModePass = false;

	//data

	$scope.schemas = [];
	Clotho.query({"schema": "Schema"}).then(function (data) {
		$scope.schemas = data;
	});

	$scope.editableList = [];
	Clotho.query({}).then(function (data) {
		var editableList = [];
		for (var i = 0; i < data.length; i++) {
			if (data[i].schema != "BuiltInSchema") {
				editableList.push(data[i]);
			}
		}
		$scope.editableList = editableList;
	});

	//functionality

	$scope.createNewObject = function () {
		Clotho.create({schema: $scope.selected}).then(function (id) {
			console.log(id);
			$scope.editable = id;
			$scope.editModePass = true;
		});
		$scope.selected = undefined;
	};

	$scope.logSelected = function () {
		console.log($scope.selected);
	};

	$scope.createNewSchema = function () {
		$scope.editModePass = true;
		$scope.editable = {schema: "ClothoSchema", language: "JSONSCHEMA"};
	};
});

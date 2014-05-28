'use strict';

angular.module('clotho.webapp')
  .controller('TestConstructionCtrl', function ($scope, $http, ConstructionSimulator) {


		/*
		$http.get('models/construction/construction_pHA581.json')
		.success(function(data) {
			$scope.demoConstruction = data;
			ConstructionSimulator.process(data).then(function (processedFile) {
				$scope.processed = processedFile;
			})
		});
		*/

		$http.get('models/construction/construction_gfp.json')
			.success(function(data) {
				$scope.demoConstruction = data;
				ConstructionSimulator.process(data).then(function (processedFile) {
					$scope.processed = processedFile;
				})
			});

		$http.get('models/construction/construction_parsed_kan.json').success(function(data) {
			$scope.parsed = data;
		});

	});

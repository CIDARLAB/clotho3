'use strict';

angular.module('clotho.webapp')
  .controller('TestConstructionCtrl', function ($scope, $http, ConstructionSimulator) {

		var files = ['pHA581', 'pSB1A2', 'basic', 'gfp'];

	  $http.get('models/construction/construction_'+files[1]+'.json')
		 .success(function(data) {
			$scope.demoConstruction = data;
			$scope.demoParse = angular.copy(data);
		});

		$http.get('models/construction/construction_'+files[3]+'.txt')
		.success(function(data) {
			$scope.demoWetLab = data;
		});


	});

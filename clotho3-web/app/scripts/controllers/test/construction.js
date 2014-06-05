'use strict';

angular.module('clotho.webapp')
  .controller('TestConstructionCtrl', function ($scope, $http, ConstructionSimulator) {

		//$http.get('models/construction/construction_pHA581.json')
	  $http.get('models/construction/construction_pSB1A2.json')
	  //$http.get('models/construction/construction_basic.json')
	  //$http.get('models/construction/construction_gfp.json')
		 .success(function(data) {
			$scope.demoConstruction = data;
			$scope.demoParse = angular.copy(data);
		});
	});

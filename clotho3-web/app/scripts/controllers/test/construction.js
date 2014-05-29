'use strict';

angular.module('clotho.webapp')
  .controller('TestConstructionCtrl', function ($scope, $http, ConstructionSimulator) {

		//$http.get('models/construction/construction_pHA581.json')
		//$http.get('models/construction/construction_parsed_kan.json')
	  $http.get('models/construction/construction_demo2.json')
	  //$http.get('models/construction/construction_basic.json')
		 .success(function(data) {
			$scope.demoConstruction = data;
		});
	});

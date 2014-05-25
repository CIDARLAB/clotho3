'use strict';

angular.module('clotho.webapp')
  .controller('TestConstructionCtrl', function ($scope, $http) {

		$http.get('models/construction/construction_parsed_kan.json').success(function(data) {
			$scope.parsed = data;
		});

	});

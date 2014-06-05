'use strict';

$clotho.extensions.controller('constructionFiles_runningCFCtrl', function($scope, $http) {

	$http.get('models/construction/construction_pSB1A2.json').success(function(data) {
		$scope.parsed = data;
	});

});
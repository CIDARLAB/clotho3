'use strict';

$clotho.extensions.controller('constructionFiles_constructionInClothoCtrl', function($scope, $timeout, $http, Clotho) {

	$http.get('models/construction/construction_pSB1A2.txt', {cache : true}).success(function (data) {
		$scope.wetlab = data;
	});

	$http.get('models/construction/construction_pSB1A2.json').success(function(data) {
		$scope.parsed = data;
	});

});
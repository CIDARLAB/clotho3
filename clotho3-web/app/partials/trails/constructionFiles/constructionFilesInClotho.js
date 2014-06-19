'use strict';

$clotho.extensions.controller('constructionFiles_constructionInClothoCtrl', function($scope, $http, Clotho) {

	$http.get('models/construction/construction_gfp.json').success(function(data) {
		$scope.parsed = data;
	});

	Clotho.get('clotho.demo.sequence.pPROBE-GFP').then(function (result) {
		$scope.demoseq = result;
	});

});
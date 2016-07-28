'use strict';

$clotho.extensions.controller('constructionFiles_constructionInClothoCtrl', function($scope, $http, Clotho) {

	$http.get('models/construction/construction_gfp.json').success(function(data) {
		$scope.parsed = data;
	});

	Clotho.get('clotho.plasmid.pPROBE-GFP[LVA]').then(function (result) {
		$scope.demoseq = result;
	});

});
'use strict';

$clotho.extensions.controller('constructionFiles_parsingCFCtrl', function($scope, $http) {

	$http.get('models/construction/construction_gfp.txt')
		.success(function(data) {
			$scope.demoWetLab = data;
		});

});
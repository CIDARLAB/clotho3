'use strict';

$clotho.extensions.controller('constructionFiles_wetLabSyntaxCtrl', function ($scope, $http) {

	$http.get('models/construction/construction_gfp.txt', {cache : true})
	.success(function (data) {
		$scope.gfpFile = data;
	});

});
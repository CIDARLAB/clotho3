'use strict';

$clotho.extensions.controller('constructionFiles_wetLabExamplesCtrl', function($scope, $http, Clotho) {

	//todo - add more versions

	$scope.files = {
		"gfp" : {},
		"pHA581" : {},
		"pSB1A2-Bca9128" : {},
		"vio" : {}
	};

	angular.forEach($scope.files, function (obj, name) {
		$http.get('models/construction/construction_' + name + '.txt', {cache : true})
			.success(function (data) {
				angular.extend(obj, {
					name : name,
					data : data
				});
			})
	});

	$scope.showFile = function (file) {
		$scope.showing = file;
	};

	$scope.showFile($scope.files.gfp);

});
'use strict';

$clotho.extensions.controller('constructionFiles_runningCFCtrl', function($scope, $http) {

	//todo - add more versions

	$scope.files = {
                "pSB1A2-Bca9128" : {},
                "vio" : {},
		"pHA581" : {}

	};

	angular.forEach($scope.files, function (obj, name) {
		$http.get('models/construction/construction_' + name + '.json', {cache : true})
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

	$scope.showFile($scope.files.gfp)

});
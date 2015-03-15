'use strict';

$clotho.extensions.controller('constructionFiles_goingThroughCtrl', function($scope, $timeout, $http, Clotho) {

	$http.get('models/construction/construction_gfp.json').then(function(data) {
		$scope.constructionFile = data.data
	});

});
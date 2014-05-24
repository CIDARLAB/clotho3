'use strict';

$clotho.extensions.controller('constructionFiles_constructionInClothoCtrl', function($scope, $timeout, $http, Clotho) {

	$http.get('models/construction/construction_parsed_kan.json').success(function(data) {
		$scope.parsed = data;
	});

});
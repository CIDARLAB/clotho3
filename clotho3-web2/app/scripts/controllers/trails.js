angular.module('clotho.webapp').controller('TrailsCtrl', function ($scope, Clotho) {
	$scope.trails = [];

	Clotho.query({schema : "Trail"}).then(function(result) {
		$scope.trails = result;
	});

	$scope.base64icon = base64icon;
});

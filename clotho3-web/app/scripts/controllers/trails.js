angular.module('clotho.webapp').controller('TrailsCtrl', function ($scope, Clotho) {
	$scope.trails = [];

	Clotho.query({schema : "org.clothocad.model.Trail"}).then(function(result) {
		$scope.trails = result;
	});
});

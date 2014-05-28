angular.module('clotho.webapp').controller('TrailsCtrl', function ($scope, $location, Clotho) {

	$scope.headerMockTrail = {
		name : "Trail Browser",
		description : "Learn Synthetic Biology with Clotho"
	};

	$scope.trails = [];
	Clotho.query({schema : "org.clothocad.model.Trail"}).then(function(result) {
		$scope.trails = result;
	});

	$scope.defaultTrailIcon = 'images/trails/trails_logo.png';

	$scope.startTrail = function (id) {
		console.log(id);
		Clotho.startTrail(id);
	};

	$scope.startTrailPage = function (id, pos) {
		$location.search('position', pos);
		Clotho.startTrail(id);
	};

	$scope.highlight = function (trail, evt) {
		//highlight is trail in object above, selected is the actual trail
		$scope.highlighted = trail;
		$scope.loading = true;
		Clotho.get(trail.id, {mute : true}).then(function (result) {
			$scope.loading = false;
			$scope.selected = result;
		});
	}
});

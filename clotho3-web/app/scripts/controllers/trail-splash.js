'use strict';

angular.module('clotho.trails')
	.controller('TrailSplashCtrl', function ($scope, $location, $http, Clotho) {

		$http.get('models/trail-splash.json', {cache : true}).success(function (data) {
			$scope.topics = data;
		});

		$scope.startTrail = function (id) {
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

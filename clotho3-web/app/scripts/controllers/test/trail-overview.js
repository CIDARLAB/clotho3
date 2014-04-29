'use strict';

angular.module('clotho.trails')
	.controller('TestTrailOverviewCtrl', function ($scope, Clotho) {
		Clotho.get('bb99191e810c19729de860fe').then(function (result) {
			$scope.trail = result;
		});

		$scope.defaultTrailIcon = 'images/trails_splash.jpg';
	});
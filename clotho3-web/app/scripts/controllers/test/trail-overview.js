'use strict';

angular.module('clotho.trails')
	.controller('TestTrailOverviewCtrl', function ($scope, $http, Clotho) {
		//Clotho.get('bb99191e810c19729de860fe').then(function (result) {
		$http.get('models/org.clothocad.trails.LearningClotho.json').success(function (data, headers) {
			$scope.trail = data;
		});
	});
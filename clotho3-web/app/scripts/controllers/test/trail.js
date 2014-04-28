'use strict';

angular.module('clotho.trails')
  .controller('TestTrailCtrl', function ($scope, Clotho) {
		Clotho.get('bb99191e810c19729de860fe').then(function (result) {
			$scope.trail = result;
		});
  });

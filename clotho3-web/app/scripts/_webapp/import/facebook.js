'use strict';

angular.module('clotho.webapp')
.controller('FacebookImportCtrl', function ($scope, Clotho, Facebook) {

	Facebook.getUser().then(function (user) {
		//got the user!
		console.log('got user');
		$scope.user = user;
	}, function () {
		//need to login
		console.log('need to login');
		$scope.showLogin = true;
	});

	$scope.login = function () {
		Facebook.login().then(function (user) {
			$scope.user = user;
		}, function (err) {
			console.log('something didn\'t work');
			//todo - we will need to re-request if denied: see the API
		});
	};

	$scope.logout = function () {
		Facebook.logout().then(function () {
			$scope.user = null;
			$scope.retrieved = null;
			$scope.createdId = null;
			$scope.showLogin = true;
		});
	};

	$scope.$watch('user.email', function (newval) {
		newval && Clotho.get(newval, {mute : true}).then(function (retrieved) {
			$scope.retrieved = retrieved;
			$scope.createdId = retrieved.id;
		});
	});

	$scope.createUser = function () {
		if ($scope.user.email) {
			var person = Facebook.convertToPersonSharable($scope.user);

			if (!$scope.retrieved) {
				Clotho.create(person).then(function (id) {
					$scope.createdId = id;
				}, function () {
					//already existed
				});
			}
		}
	}

});

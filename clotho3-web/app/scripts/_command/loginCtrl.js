angular.module('clotho.commandbar')
	.controller('loginCtrl', function ($scope, $timeout, Clotho) {

		$scope.notification = {};
		$scope.cred = {username: "", password: ""};

		$scope.login = function () {
			Clotho.login($scope.cred.username, $scope.cred.password)
			.then(function (result) {
				$scope.notification = {class: "alert-success", message: "Log in Success"};
			}, function (err) {
				$scope.notification = {class: "alert-danger", message: "Log in Error"};
				$scope.cred = {username: "", password: ""};
			});
		};

	});

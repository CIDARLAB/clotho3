angular.module('clotho.commandbar')
	.controller('loginCtrl', function ($scope, $timeout, Clotho) {

		$scope.notification = {};
		$scope.cred = {username: "", password: ""};

		$scope.login = function () {
			Clotho.login($scope.cred.username, $scope.cred.password).then(function (result) {
				console.log('run login', result);
				if (!!result) {
					$scope.notification = {class: "alert-success", message: "Log in Success"};
				} else {
					$scope.notification = {class: "alert-danger", message: "Log in Error"};
					$scope.cred = {username: "", password: ""};
				}
			});
		};

	});

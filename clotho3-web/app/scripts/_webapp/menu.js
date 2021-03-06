'use strict';

angular.module('clotho.webapp')
  .controller('MenuCtrl', function ($scope, $location, $timeout, Collector, PubSub, hotkeys, Clotho, $modal) {

		$scope.modes = [
			{"name" : "Editor", "path" : "/editor"},
			{"name" : "Trails", "path" : "/trails"},
			{"name" : "Browser", "path" : "/browser"}
		];

		$scope.$watch(function () {
			return $location.path();
		}, function (newValue, oldValue) {
			if (!$scope.modes) return;
			angular.forEach($scope.modes.items, function(mode, num) {
				var regexp = new RegExp('^' + mode.path + '.*$', ['i']);
				if (regexp.test(newValue)) {
					mode.class = "active";
				} else {
					//remove class
					mode.class = ""
				}
			});
		});

		//initial check of url, not picked up in watch
		$timeout(function() {
			var path = $location.path();
			angular.forEach($scope.modes.items, function(mode, num) {
				var regexp = new RegExp('^' + mode.path + '.*$', ['i']);
				mode.class = (regexp.test(path)) ? 'active' : '';
			});
		}, 0);


		//do hrefs or this make more sense?
		$scope.goToPage = function(mode) {
			$location.path($scope.modes[mode]);
			//note: $location.path().replace() will avoid creation of page in history
		};

		$scope.clearStorage = function() {Collector.clearStorage()};
		$scope.logListeners = function() {PubSub.logListeners()};

		$scope.loggedIn = false;
		$scope.showLogin = function() {
			$modal.login().result.then(function(result) {
				console.log(result);
				if (result) {
					$scope.username = result;
					$scope.loggedIn = true;
				}
			});
		};

		$scope.hideNavBar = true;
		$scope.showMenu = function() {
			console.log('showing');
			$scope.hideNavBar = !$scope.hideNavBar;
		};

		hotkeys.add('alt+ctrl+,', 'Show Route Menu', $scope.showMenu, true);

  });

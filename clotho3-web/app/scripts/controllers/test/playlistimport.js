'use strict';

angular.module('clotho.webapp')
	.controller('TestPlaylistimportCtrl', function ($scope, Youtube, $q) {

		$scope.playlistId = 'PL2aPXzks-TgO0k9PhT__NSh2x6HNimaOy';

		$scope.$watch('playlistId', function (newval) {
			Youtube.playlistItems(newval)
			.then(function (result) {
				 $scope.playlistInfo = result;
			});
			Youtube.playlistInfo(newval)
			.then(function (result) {
				$scope.playlistItems = result;
			});
			Youtube.playlistToTrail(newval).then(function (result) {
				$scope.playlistTrail = result;
			});
		});

		//testing api get

		$scope.$watch('search', function (val)  {
			if (val) {
				Youtube.videoSearch(val).then(function (data) {
					$scope.searchResult = data;
				});
			} else {
				$scope.searchResult = '';
			}
		});

	});

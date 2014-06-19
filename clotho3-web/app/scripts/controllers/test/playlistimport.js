'use strict';

angular.module('clotho.webapp')
	.controller('TestPlaylistimportCtrl', function ($scope, Youtube, Clotho, $q) {

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

		$scope.create = function () {
			if (!angular.isEmpty($scope.playlistTrail)) {
				Clotho.create($scope.playlistTrail);
			}
		};

		/*
		These have already been imported. This is just to illustrate flow at this point

		$scope.playlists = [
			'PL2aPXzks-TgMbE-b15ezcHQmGVTfMuHcZ',
			'PL2aPXzks-TgO0k9PhT__NSh2x6HNimaOy',
			'PL2aPXzks-TgNkKi2lUWL65MpPsRBCb5q5',
			'PL2aPXzks-TgPX7HU4pogNUSBN9OK2PhBK',
			'PL2aPXzks-TgOtngzKdocz4P50TCncRyPy',
			'PL2aPXzks-TgMvy4PG5wwbrw9VKBuowOR8',
			'PL2aPXzks-TgNP1FuAv0v-fS6UhCjy1waN',
			'PL2aPXzks-TgN_Ztp-bRkzbtbi5Ncsh71t',
			'PL2aPXzks-TgMo_x6XhYmfoZDezm4m-DcL',
			'PL2aPXzks-TgOj55WalLbDUTfmuhZvYVmW',
			'PL2aPXzks-TgNF957XOSBB3rPdFbtvsWYI',
			'PL2aPXzks-TgOOPfxneYJ3pBwTlqzv9Ezn',
			'PL2aPXzks-TgNlrgM8BK3jD_aiyjoXOsRI'
		];

		$scope.createAll = function () {
			angular.forEach($scope.playlists, function (id) {
				Youtube.playlistToTrail(id).then(function (result) {
					Clotho.create(result);
				});
			})
		};
		*/

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

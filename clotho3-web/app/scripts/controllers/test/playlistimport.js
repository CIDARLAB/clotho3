'use strict';

angular.module('clotho.webapp')
	.controller('TestPlaylistimportCtrl', function ($scope, Youtube, $q) {

		$scope.playlistId = 'PL2aPXzks-TgO0k9PhT__NSh2x6HNimaOy';

		$scope.$watch('playlistId', function (newval) {
			$q.all({
				playlistItems : Youtube.playlistItems(newval),
				playlistInfo : Youtube.playlistInfo(newval)
			}).then(function (resultObj) {
					$scope.showerror = false;
					$scope.playlistInfo = resultObj.playlistInfo;
					$scope.playlistItems = resultObj.playlistItems;
					$scope.result = parsePlaylist($scope.playlistInfo, $scope.playlistItems);
				},
				//error handling
				function (err) {
					console.log('error');
					$scope.showerror = true;
					$scope.playlist = {};
					$scope.result = {};
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

		function parsePlaylist (info, items) {

			//parse out info we want

			var result = {
				name : info.snippet.title,
				description : info.snippet.description,
				created : new Date(info.snippet.publishedAt).valueOf(),
				icon : info.snippet.thumbnails.standard.url
			};

			//parse out videos

			result.videos = [];
			_.each(items, function (item) {
				result.videos.push({
					videoId : item.contentDetails.videoId,
					name : item.snippet.title,
					description : item.snippet.description,
					icon : item.snippet.thumbnails.standard.url
				})
			});

			return result;
		}

	});

angular.module('clotho.webapp')
.controller('ImportCtrl', function ($scope, $location, Clotho) {

		//things we know how to import
		$scope.importable = [
			{
				name : "Youtube Playlist",
				route : '/import/youtubePlaylist'
			},
			{
				name : "ApE Files",
				route : '/import/ape'
			}
		];

		$scope.goToImportable = function (importable) {
			$location.path(importable.route);
		}
});

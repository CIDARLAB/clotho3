angular.module('clotho.trails')
	.service('Youtube', function($http, $rootScope, $q, $timeout, $window) {

		//load the API and return it as promise resolution
		var api_ready = $q.defer();
		$window.onYouTubeIframeAPIReady = function () {
			api_ready.resolve(YT);
		};
		$clotho.extensions.script('//www.youtube.com/iframe_api');

	/**
	 * @description Given a URL (youtube.com, youtu.be, watch, embed, etc.), extracts the youtube VideoID. Passing in a VideoId will work.
	 * @source http://stackoverflow.com/a/10315969/624466
	 * @note assumes length of 11 characters. Youtube may change this in the future.
	 *
	 * @param {string} url URL for the video, containing the ID
	 * @returns {string} videoId
	 */
	var extract = function extractYoutubeID (url) {
		var regex = /^(?:https?:\/\/)?(?:www\.)?(?:youtu\.be\/|youtube\.com\/(?:embed\/|v\/|watch\?v=|watch\?.+&v=))((\w|-){11})(?:\S+)?$/;
		return (url.match(regex) || url.match(/((\w|-){11})/)) ? RegExp.$1 : false;
	};

	/**
	 * @description Given a videoId, retrieves the youtube feeds API information about it
	 *
	 * @param {string} videoId ID of the video
	 * @param {string} size Thumbnail size, youtube names (default, hqdefault, mqdefault, 0, 1, 2, 3). default is 'default' (120x90)
	 * @returns {Promise} thumbnail URL
	 */
	var thumbnail = function generateYoutubeThumbnailUrl(videoId, size) {
		size = size || "default";

		return "https://img.youtube.com/vi/"+ videoId + "/" + size + ".jpg";
	};

	/**
	 * @description Given a valid videoId, retrieves the youtube feeds API information about it
	 *
	 * @param {string} videoId
	 * @returns {Promise} data (JSON) retrieved from youtube feeds API
	 */
	var info = function getYoutubeInfo (videoId) {
		return $http.get('https://gdata.youtube.com/feeds/api/videos/'+videoId+'?&alt=json')
			.then(function (data) {
				return data.data
			})
	};

	/**
	 * @description Given a valid playlist, retrieves the youtube feeds API information about it
	 *
	 * @param {string} playlistId
	 * @returns {Promise} data (JSON) retrieved from youtube feeds API
	 */
	var playlist = function getYoutubePlaylist (playlistId) {
		return $http.get('http://gdata.youtube.com/feeds/api/playlists/'+playlistId+'?alt=json')
			.then(function (data) {
				return data.data
			})
	};

	return {
		apiPromise : api_ready.promise,
		extract : extract,
		thumbnail : thumbnail,
		info : info,
		playlist : playlist
	}
});
angular.module('clotho.trails')
	.config(function() {
		//todo - move to $script or $clotho.extensions.script()
		//load the Youtube API
		var tag = document.createElement('script');
		tag.src = "//www.youtube.com/iframe_api";
		var firstScriptTag = document.getElementsByTagName('script')[0];
		firstScriptTag.parentNode.insertBefore(tag, firstScriptTag);
	})
	.service('Youtube', function($http, $rootScope, $q, $timeout) {

		var api_ready = $q.defer();
		$rootScope.$watch(function() {
			return YT.loaded == 1;
		},function(newval) {
			if (!!newval) {
				api_ready.resolve();
			} else {
				//note - force a digest after some time if api not ready, hopefully ready by then
				$timeout(angular.noop, 500);
			}
		});

	/**
	 * @description Given a URL (youtube.com, youtu.be, watch, embed, etc.), extracts the youtube VideoID. Passing in a VideoId will work.
	 * @source http://stackoverflow.com/a/10315969/624466
	 * @note assumes length of 11 characters. Youtube may change this in the future.
	 *
	 * @param {string} url
	 * @returns {string} videoId
	 */
	var extract = function extractYoutubeID (url) {
		var regex = /^(?:https?:\/\/)?(?:www\.)?(?:youtu\.be\/|youtube\.com\/(?:embed\/|v\/|watch\?v=|watch\?.+&v=))((\w|-){11})(?:\S+)?$/;
		return (url.match(regex) || url.match(/((\w|-){11})/)) ? RegExp.$1 : false;
	};

	//can use youtube names (default, hqdefault, mqdefault, 0, 1, 2, 3)
	//defaults to default (120 x 90)
	var thumbnail = function generateYoutubeThumbnailUrl(videoId, size) {
		size = size || "default";

		return "https://img.youtube.com/vi/"+ videoId + "/" + size + ".jpg";
	};

	var info = function getYoutubeInfo (videoId) {
		return $http.get('https://gdata.youtube.com/feeds/api/videos/'+videoId+'?v=2&prettyprint=true&alt=jsonc')
			.then(function (data) {
				return data.data
			})
	};

	return {
		readyPromise : api_ready.promise,
		extract : extract,
		thumbnail : thumbnail,
		info : info
	}
});
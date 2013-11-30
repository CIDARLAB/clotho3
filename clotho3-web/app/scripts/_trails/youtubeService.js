angular.module('clotho.trails').service('YoutubeService', function($http) {

	/**
	 * @description Given a URL (youtube.com, youtu.be, watch, embed, etc.), extracts the youtube VideoID. Passing in a VideoId will work.
	 * @source http://stackoverflow.com/a/10315969/624466
	 * @note assumes length of 11 characters. Youtube may change this in the future.
	 *
	 * @param {string} url
	 * @returns {string} videoId
	 */
	var extract_youtube = function(url) {
		var regex = /^(?:https?:\/\/)?(?:www\.)?(?:youtu\.be\/|youtube\.com\/(?:embed\/|v\/|watch\?v=|watch\?.+&v=))((\w|-){11})(?:\S+)?$/;
		return (url.match(regex) || url.match(/((\w|-){11})/)) ? RegExp.$1 : false;
	};

	//can use youtube names (default, hqdefault, mqdefault, 0, 1, 2, 3)
	//defaults to default (120 x 90)
	var youtubeThumbnail = function(videoId, size) {
		size = size || "default";

		return "https://img.youtube.com/vi/"+ videoId + "/" + size + ".jpg";
	};

	var youtubeInfo = function(videoId) {
		return $http.get('https://gdata.youtube.com/feeds/api/videos/'+videoId+'?v=2&prettyprint=true&alt=jsonc')
			.then(function (data) {
				return data.data
			})
	};

	return {
		extract_youtube : extract_youtube,
		youtubeThumbnail : youtubeThumbnail,
		youtubeInfo : youtubeInfo
	}
});
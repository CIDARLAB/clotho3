angular.module('clotho.youtube')
	.service('Youtube', function($http, $rootScope, $q, $timeout, $window) {

		//todo - system to manage keys and change depending on build, without exposing testing keys
		//other keys at https://console.developers.google.com/project/apps~clotho-bio/apiui/credential
		//key only works at *.synbiotrails.org/*

		var PUBLIC_KEY = 'AIzaSyDcv-f_kDYSvRMpsfKn9gaWo3Njzr9QmDU';
		var url_base = 'https://www.googleapis.com/youtube/v3';

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
	var videoThumbnail = function generateYoutubeThumbnailUrl(videoId, size) {
		size = size || "default";

		return "https://img.youtube.com/vi/"+ videoId + "/" + size + ".jpg";
	};

	var videoSearch = function youtubeSearch (term) {
		return $http.get(url_base + '/search', {
			params: {
				key: PUBLIC_KEY,
				type: 'video',
				maxResults: '15',
				part: 'id, snippet',
				fields: 'items/id,items/snippet/title,items/snippet/description,items/snippet/thumbnails/default,items/snippet/channelTitle',
				q: term
			}
		}).then(function (data) {
			return data;
		});
	};

	/**
	 * @description Given a valid videoId, retrieves the youtube feeds API information about it
	 *
	 * @param {string} videoId
	 * @returns {Promise} data (JSON) retrieved from youtube feeds API
	 */
	var videoInfo = function getYoutubeInfo (videoId) {
		return $http.get('https://gdata.youtube.com/feeds/api/videos/'+videoId+'?&alt=json')
			.then(function (data) {
				return data.data.entry
			})
	};

	/**
	 * @description
	 * Given a valid playlist id, retrieves information about playlist.
	 *
	 * We parse this down to the single list item, so you get a [youtube playlist](https://developers.google.com/youtube/v3/docs/playlists)
	 *
	 * Doesn't retrieve contents - see playlistItems for that.
	 *
	 * @param {string} playlistId
	 * @returns {Promise} data (JSON) retrieved from youtube feeds API
	 */
	var playlistInfo = function getYoutubePlaylistInfo (playlistId) {
		return $http.get('https://www.googleapis.com/youtube/v3/playlists', {
			params: {
				key: PUBLIC_KEY,
				part: 'id, snippet, status, contentDetails',
				id: playlistId
			},
			cache : true
		})
		.then(function (data) {
			return data.data.items[0];
		})
	};

	/**
	 * @description
	 * Given a valid playlist, retreives information about contents.
	 *
	 * Note that this will not handle playlists with more than 50 items
	 *
	 * https://developers.google.com/youtube/v3/docs/playlistItems
	 *
	 * @param {string} playlistId
	 * @returns {Promise} data (JSON) retrieved from youtube feeds API
	 */
	var playlistItems = function getYoutubePlaylist (playlistId) {
		return $http.get(url_base + '/playlistItems', {
			params: {
				key: PUBLIC_KEY,
				maxResults: '50',
				part: 'id, snippet, status, contentDetails',
				playlistId: playlistId
			},
			cache: true
		})
		.then(function (data) {
			return data.data.items;
		})
	};

	var playlistToTrail = function (playlistId) {
		return $q.all({
			playlistItems : playlistItems(playlistId),
			playlistInfo : playlistInfo(playlistId)
		})
		.then(function (resultObj) {
			var info = resultObj.playlistInfo;
			var items = resultObj.playlistItems;

			var result = {
				name : info.snippet.title,
				schema : "org.clothocad.model.Trail",
				id : "org.clothocad.trails.youtube." + (info.snippet.title).replace(/\s+/g, ''),
				about : {
					help : info.snippet.title + " was retrieved from Youtube from the channel " + info.snippet.channelTitle + ", and was created on " + (new Date(info.snippet.publishedAt)).toLocaleDateString() + ".",
					contents : [{
						type : "text",
						params: info.snippet.description
					}]
				},
				created : new Date(info.snippet.publishedAt).valueOf(),
				youtubePlaylistId: playlistId,
				icon : !!info.snippet.thumbnails.standard ?
								 info.snippet.thumbnails.standard.url :
								 info.snippet.thumbnails.high.url,
				contents: [
					{
						chapter: 'Youtube Playlist',
						description : 'This section was imported from youtube',
						pages: []
					}
				]
			};

			//parse out videos
			_.each(items, function (item) {

				var newPage = {
					title : item.snippet.title,
					icon : "video",
					contents : []
				};

				if (item.snippet.description && item.snippet.description.length) {
					newPage.contents.push({
						type: "text",
						params : item.snippet.description
					});
				}

				newPage.contents.push({
					type: "video",
					params: {
						id: item.contentDetails.videoId
					}
				});

				result.contents[0].pages.push(newPage)
			});

			return result;
		});
	};

	return {
		apiPromise : api_ready.promise,
		extract : extract,

		videoSearch: videoSearch,
		videoThumbnail : videoThumbnail,
		videoInfo : videoInfo,

		playlistInfo : playlistInfo,
		playlistItems : playlistItems,

		playlistToTrail : playlistToTrail
	}
});
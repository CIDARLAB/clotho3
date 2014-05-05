/**
 * @name youtube
 * @ngdoc directive
 *
 * @params
 * videoId {string}
 * params {object} Youtube video parameters for construction, as defined at https://developers.google.com/youtube/iframe_api_reference. Note onStateChange will be overwritten, but other events can be passed here
 * autoplay {boolean}
 * startMini {boolean}
 * onPlay {function}
 * onPause {function}
 * onComplete {function}
 */

angular.module('clotho.trails')
	.directive('youtube', function (Trails, Youtube, $compile, $timeout, $http) {

		return {
			restrict: 'EA',
			replace: true,
			scope: {
				videoId: '@youtube',
				params: '=?',
				autoplay: '@?', //can also define on params directly, as 'autoplay'
				startMini: '@?', //can also define on params directly, as 'mini'
				onComplete: '&?',
				onPlay: '&?',
				onPause: '&?'
			},
			link: function youtubeDirectiveLink(scope, element, attrs) {

				//defaults
				var defaults = {
					width: 700,
					height: 525,
					border: 0,
					autoplay: false,
					mini: false,
					videoId: scope.videoId,
					playerVars: {
						autoplay: (!scope.autoplay && !scope.startMini) ? 0 : 1,
						autohide: 1,
						rel: 0
					},
					events: {}
				};
				scope.params = angular.extend(defaults, scope.params);

				if (!scope.params.videoId) return;

				//pull out of params so don't need to declare as attribute
				scope.autoplay = angular.isDefined(scope.autoplay) ?
					scope.autoplay : scope.params.autoplay;
				scope.startMini = angular.isDefined(scope.startMini) ?
					scope.startMini : scope.params.mini;

				scope.params.events.onStateChange = function (event) {
					if (event.data == 0) {
						scope.$apply(scope.onComplete());
					}
					else if (event.data == 1) {
						scope.$apply(scope.onPlay());
					}
					else if (event.data == 2) {
						scope.$apply(scope.onPause());
					}
				};

				scope.convertToPlayer = function () {
					createYoutubePlayer();
				};

				//init()
				if (!!scope.startMini && scope.startMini != 'false') {
					scope.miniThumb = Youtube.thumbnail(scope.videoId, 'mqdefault');

					Youtube.info(scope.videoId).then(function (json) {
						scope.miniInfo = json.data;
						scope.miniInfo.durationFormatted = (Math.floor(scope.miniInfo.duration / 60) + ":" + ((scope.miniInfo.duration % 60 < 10) ? '0' : '') + (scope.miniInfo.duration % 60));
					});

					$http.get('views/_trails/youtubeThumbnail.html')
						.success(function (data, headers) {
							element.html($compile(data)(scope));
						})
						.error(function (data, headers) {
							//this shouldn't happen
						});

				} else {
					element.html('<div class="youtubeLoading" width="' + scope.params.width + '" height="' + scope.params.height + '"></div>');
					createYoutubePlayer()
				}

				function createYoutubePlayer() {
					Youtube.apiPromise.then(function (YT) {
						scope.player = new YT.Player(element[0], scope.params);
						angular.element(scope.player.a).addClass('youtubePlayer');
					});
				}

				scope.$watch('videoId', function (newval, oldval) {
					if (newval != oldval) {
						scope.params.videoId = newval;
						createYoutubePlayer()
					}
				});
			}
		}
	});
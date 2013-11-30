angular.module('clotho.trails').directive('youtube', function(Trails, YoutubeService, $compile, $timeout) {
	//note - requires youtube iFrame API be present
	//todo - rewrite so youtube API loaded in this directive, match protocol (http / https)

	return {
		restrict : 'EA',
		replace: true,
		scope: {
			videoId : '@youtube',
			params : '=?',
			autoplay: '@?', //can also define on params directly, as 'autoplay'
			startMini: '@?', //can also define on params directly, as 'mini'
			onComplete : '&?'
		},
		compile: function compile(tElement, tAttrs, transclude) {
			return {
				pre: function preLink(scope, element, attrs) {

					//console.log(element);

				},
				post: function postLink(scope, element, attrs) {

					if (!scope.videoId) return;

					var videoInfo = YoutubeService.youtubeInfo(scope.videoId);

					//todo - extract dimensions from youtube info
					//todo - CSS to center if not 700px wide

					//defaults
					var defaults = {
						width : 700,
						height : 525,
						border: 0,
						autoplay: false,
						mini: false,
						videoId : scope.videoId,
						playerVars : {
							autoplay : (!scope.autoplay && !scope.startMini) ? 0 : 1,
							autohide : 1,
							rel : 0
						},
						events : {}
					};
					scope.params = angular.extend(defaults, scope.params);

					//pull out of params so don't need to declare as attribute
					scope.autoplay = angular.isDefined(scope.autoplay) ?
						scope.autoplay : scope.params.autoplay;
					scope.startMini = angular.isDefined(scope.startMini) ?
						scope.startMini : scope.params.mini;

					scope.params.events.onStateChange = function (event) {
						if (event.data == 0) {
							scope.$apply(scope.onComplete());
						}
					};

					scope.convertToPlayer = function() {
						createYoutubePlayer();
					};

					if (!!scope.startMini && scope.startMini != 'false') {
						scope.miniThumb = YoutubeService.youtubeThumbnail(scope.videoId, 'mqdefault');

						videoInfo.then(function(json) {
							scope.miniInfo = json.data;
							scope.miniInfo.durationFormatted = (Math.floor(scope.miniInfo.duration/60) + ":" + ((scope.miniInfo.duration%60 < 10) ? '0' : '') + (scope.miniInfo.duration%60));
						});

						var thumbnailHTML = '<div class="row-fluid" style="margin-bottom: 15px">' +
							'<div class="thumbnail clearfix">' +
							'<img class="span5" ng-src="{{miniThumb}}">' +
							'<div class="span7 caption">' +
							'<h5 style="margin-top:5px">{{ miniInfo.title }}</h5>' +
							'<p style="overflow:hidden; display: -webkit-box; -webkit-line-clamp: 3; -webkit-box-orient: vertical; max-height: 4.5em">{{ miniInfo.description | limitTo:300 }}</p>' +
							'<a class="btn btn-primary" ng-click="convertToPlayer()">Watch Video {{ "(" + miniInfo.durationFormatted +")" }}</a>' +
							'</div>' +
							'</div>' +
							'</div>';

						element.html($compile(thumbnailHTML)(scope));
					} else {
						//todo - convert to being a class
						element.html('<img src="../../images/assets/ajax-loader.gif" />');
						InitialLoadCreateYoutubePlayer()
					}

					function InitialLoadCreateYoutubePlayer () {
						if (YT.loaded == 1) {
							createYoutubePlayer();
						}
						else {
							scope.$watch(function() {
								return YT.loaded == 1
							}, function(newval, oldval) {
								if (!!newval) {
									console.log('youtube player API ready - setting video');
									createYoutubePlayer();
								}
							});

							//future - once load youtube API in YoutubeService, remove this hack
							$timeout(function() {
								console.log('hopefully youtube is loaded now... otherwise above code will run on next $digest');
							}, 500);

						}
					}

					function createYoutubePlayer() {
						scope.player = new YT.Player(element[0], scope.params);
						//HACK
						$(scope.player.a).css({'box-sizing' : 'border-box', 'box-shadow' : '0px 0px 8px rgba(0, 0, 0, 0.25)'})
					}
				}
			}
		}
	}
});
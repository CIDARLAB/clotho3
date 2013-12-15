angular.module('clotho.trails').directive('youtube', function(Trails, Youtube, $compile, $timeout) {
	//note - requires youtube iFrame API be present - loaded in Youtube service

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

				},
				post: function postLink(scope, element, attrs) {

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
					};

					scope.convertToPlayer = function() {
						createYoutubePlayer();
					};

					if (!!scope.startMini && scope.startMini != 'false') {
						scope.miniThumb = Youtube.thumbnail(scope.videoId, 'mqdefault');

						Youtube.info(scope.videoId).then(function(json) {
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
						//todo - move to class (loading class)
						element.html('<img src="../../images/assets/ajax-loader.gif" />');
						InitialLoadCreateYoutubePlayer()
					}

					function InitialLoadCreateYoutubePlayer () {
						Youtube.readyPromise.then(function() {
							createYoutubePlayer();
						});
					}

					function createYoutubePlayer() {
						scope.player = new YT.Player(element[0], scope.params);
						$(scope.player.a).addClass('youtubePlayer');
					}

					scope.$watch('videoId', function(newval, oldval) {
						if (newval != oldval) {
							scope.params = newval;
							createYoutubePlayer()
						}
					})
				}
			}
		}
	}
});
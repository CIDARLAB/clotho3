'use strict';

//todo - add back in base64icon (from orig module)

Application.Trails.service('Trails',
	['Clotho', 'YoutubeService', '$q',
	function(Clotho, YoutubeService, $q) {

	var compile = function TrailCompile(trail) {

		//If pass by reference (depending on Clotho.get() ) need to copy to don't edit in dependencies
		//trail = angular.copy(trail);

		var transcludes = trail.dependencies || null,
			deferred = $q.defer();

		if (!transcludes && !trail.mixin) {
			deferred.resolve(trail);
			return deferred.promise;
		}

		var final_contents = [],
			promises = [];

		//get the transcluded trails... will be fast if in collector already
		(transcludes) && angular.forEach(transcludes, function(id) {
			promises.push(Clotho.get(id));
		});

		var mixins = (trail.mixin) ? Application.mixin(trail.mixin) : $q.when();
		mixins
			.then(function() {
				return $q.all(promises)
			})
			//after download all, pluck out the chapters we need
			.then(function (downloads) {

				//reorganize transcludes so can reference by id
				transcludes = {};
				angular.forEach(downloads, function(transclude) {
					transcludes[transclude.id] = transclude;
				});

				//iterate through trail, pushing in chapters
				angular.forEach(trail.contents, function (chapter, ind) {

					if (typeof chapter.transclude == 'undefined') {
						final_contents.push(chapter);
					} else {
						//chapters to include :
						var chapterID = chapter.transclude.id,
							chapterNum = chapter.transclude.chapters;

						if ((chapterNum == "all") || (typeof chapterNum == 'undefined')) {
							for (var i = 0; i < transcludes[chapterID]['contents'].length; i++) {
								final_contents.push(transcludes[chapterID]['contents'][i]);
							}
						} else {
							var startStop = chapterNum.split("-");
							if (startStop.length == 1) {
								final_contents.push(transcludes[chapterID]['contents'][startStop[0]]);
							} else {
								if (startStop[0] > startStop[1])
									return "wrong format - start must be smaller than end";

								for (var i = startStop[0]; i <= startStop[1]; i++) {
									final_contents.push(transcludes[chapterID]['contents'][i])
								}
							}
						}
					}
				});

				trail.contents = final_contents;
				deferred.resolve(trail);
			});

		return deferred.promise;
	};

    //in form <Chapter>-<Page>
    var extractPage = function(Trail, indices) {
        var pos = indices.split("-");
        var page = Trail.contents[pos[0]]['pages'][pos[1]];
        return page;
    };

    //bring back in logic from trail-module.js
    var calcNextPage = function(Trail, oldpos) {
        oldpos = (typeof oldpos != 'undefined') ? oldpos.split("-") : [0, -1];
        var newpos;

        if (typeof Trail.contents[oldpos[0]]['pages'][+oldpos[1] + 1] != 'undefined')
            newpos = oldpos[0] + '-' + (+oldpos[1] + 1);
        else if (typeof Trail.contents[+oldpos[0] + 1]['pages'] != 'undefined')
            newpos = (+oldpos[0] + 1) + '-' + 0;
        else {
            return;
        }
        return newpos;
    };

    var calcPrevPage = function(Trail, oldpos) {
        if (oldpos == '0-0') return;

        oldpos = (typeof oldpos != 'undefined') ? oldpos.split("-") : [0, 1];
        var newpos;

        if (typeof Trail.contents[oldpos[0]]['pages'][+oldpos[1] - 1] != 'undefined')
            newpos = oldpos[0] + '-' + (+oldpos[1] - 1);
        else if (typeof Trail.contents[+oldpos[0] - 1]['pages'])
            newpos = (+oldpos[0] - 1) + '-' + (Trail.contents[+oldpos[0] - 1]['pages'].length - 1);
        else {
            return;
        }
        return newpos;
    };



	var favorite = function(id) {
		console.log("favorite trail with id: " + id);
	};

	// PERSISTENCE

	var persist = {};

	return {
		compile : compile,
        extractPage : extractPage,
        calcNextPage : calcNextPage,
        calcPrevPage : calcPrevPage,
		share : Clotho.share,
		favorite : favorite,
		persist : persist
	}
}]);



Application.Trails.controller('TrailMainCtrl', ['$scope', 'Clotho', function($scope, Clotho) {

	$scope.trails = [];

	Clotho.query({schema : "Trail"}).then(function(result) {
		$scope.trails = result;
	});

	$scope.base64icon = base64icon;
}]);

Application.Trails.controller('TrailDetailCtrl', ['$scope', '$route', 'Clotho', 'Trails', '$http', '$timeout', '$templateCache', '$compile', '$keypress', '$q', '$controller', '$window', function($scope, $route, Clotho, Trails, $http, $timeout, $templateCache, $compile, $keypress, $q, $controller, $window) {

	//inherited from $routeProvider.resolve clause in application.js
	$scope.id = $route.current.params.id;
	$scope.trail = $route.current.locals.trail;

	//kickoff
	//$scope.content = $scope.trail.description
	$timeout(function() { $scope.activate('0-0') });


	$scope.activate = function(indices) {
        //if passed nothing
        if (!indices || !angular.isString(indices)) return;

		//don't activate already active one
		if ($scope.current == indices) return;
		$scope.current = indices;

        $scope.currentPage = Trails.extractPage($scope.current);
		//$scope.content = Trails.loadPage(page, $scope);
	};

	$scope.home = function() {
		$scope.content = $scope.trail.description;
		$scope.current = undefined;
	};

	$scope.favorite = function() {
		//future - better checking, do initial check
		$scope.favorited = !$scope.favorited;
		Trails.favorite($scope.id);
	};

	$scope.share = function() {
		Trails.share($scope.id)
	};

	$scope.next = function() {
		$scope.activate(Trails.calcNextPage($scope.trail, $scope.current));
	};

	$scope.prev = function() {
		$scope.activate(Trails.calcPrevPage($scope.trail, $scope.current));
	};

	$keypress.on('keydown', {'alt-right' : 'next()', 'alt-left' : 'prev()'}, $scope);

	$scope.base64icon = base64icon;

}]);

Application.Trails.directive('trailQuiz', ['$http', '$templateCache', '$compile', 'Clotho', '$interpolate', '$q', function($http, $templateCache, $compile, Clotho, $interpolate, $q) {
	return {
		restrict: "EA",
		require: 'ngModel',
		scope: {
			quiz: '=ngModel',
			gradeCallback : '=?',
			advance : '&?'
		},

		compile: function compile(tElement, tAttrs, transclude) {
			return {
				pre: function preLink(scope, element, attrs) {

					//todo -- rethink extending scope with whole quiz (i.e. want dictionary, but maintain quiz namespace?)
					angular.extend(scope, scope.quiz);
					//can't use $interpolate - need to maintain bindings
					//can't compile with scope.quiz - not a scope object - and can't create isolate because bindings not maintained
					// note - see also grade and retry functions, and load.quiz
					scope.quiz.question = $compile('<h5>' + scope.quiz.question + '</h5>')(scope);

					//console.log(scope.quiz.question);

					$http.get('partials/trails/quiz/' + scope.quiz.type + '-partial.html', {cache: $templateCache})
						.success(function (data) {
							element.html($compile('<div class="quiz">' + data + '</div>')(scope));
						})
						.error(function(data, status, headers, config) {
							element.html('<p>Template could not be found...</p>' + JSON.stringify(scope.quiz));
						});

				},
				post: function postLink(scope, element, attrs) {

					scope.createEmptyAnswer = function(quiz, value) {
						value = (typeof value != 'undefined') ? value : false;
						scope.quiz.answer = new Array(quiz.options.length);
						for (var i = 0; i < scope.quiz.answer.length; i++) {
							scope.quiz.answer[i] = value;
						}
					};

					scope.answerUndefined = function(quiz) {
						return (typeof quiz.answer == 'undefined' || quiz.answer === '');
					};

					scope.submitQuestion = function(quiz) {
						Clotho.gradeQuiz(quiz.questionValue, quiz.answer, quiz.answerGenerator).then(function (data) {
							console.log('gradeQuiz result: ' + data);
							scope.quiz.submitted = true;
							scope.quiz.response = {};
							scope.quiz.response.result = data;
							//console.log(scope.gradeCallback);
							scope.gradeCallback(data);
						});
					};

					scope.resetQuiz = function () {
						scope.quiz.submitted = false;
						scope.quiz.response = null;
						scope.quiz.answer = null;
					};

					scope.retryQuiz = function () {
						if (!scope.quiz.retry) return;

						var promises = {},
							deferred = $q.defer();

						angular.forEach(scope.quiz.retry, function(value, key) {
							promises[key] = Clotho.submit(value).then(function (result) {
								return result;
							});
						});

						$q.all(promises).then(function(completed) {
							//todo - ugly -- get everything into quiz object
							angular.extend(scope.quiz, completed);
							angular.extend(scope, scope.quiz);
							scope.resetQuiz();
							console.log(scope);
							deferred.resolve();
						});

						return deferred.promise;
					}

				}
			}
		}
	}
}]);

Application.Trails.service('YoutubeService', ['$http', function($http) {

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
}]);

Application.Trails.directive('youtube', ['Trails', 'YoutubeService', '$compile', '$timeout', function(Trails, YoutubeService, $compile, $timeout) {
	//note - requires youtube iFrame API be present
	//todo - rewrite so youtube API loaded in this directive

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



					//todo - center if width < 700

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
						element.html('<img src="/images/ajax-loader.gif" />');
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
}]);
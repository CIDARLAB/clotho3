'use strict';

//todo - add back in base64icon

Application.Trails.service('Trails', ['Clotho', 'YoutubeService', '$q', '$dialog', '$http', '$timeout', '$templateCache', '$compile', '$controller', function(Clotho, YoutubeService, $q, $dialog, $http, $timeout, $templateCache, $compile, $controller) {







	// fields that are handled:
	// backend: CSS (url), mixin (array|url), script (array|url), onload (array|url), controller (name, must be mixed in)
	// content: text (html), video (object), template (url), quiz (object), markdown (text)
	// @param {TrailPage} page
	var loadPage = function TrailLoadPage(page, $scope) {

		//future in ng-1.2.x, use notify callbacks for updates
		//todo - error callbacks

		if (!!page.dictionary) {
			angular.extend($scope, page.dictionary);
		}

		return Application.css(page.css)
		.then(function() {
			return Application.mixin(page.mixin)
		})
		.then(function() {
			return Application.script(page.script)
		})
		.then(function (){
			var promises = [];
			angular.forEach(page.contents, function (component) {
				promises.push(loadPageComponent(component.type, component.params))
			});
			return promises;
		})
		.then(function(content) {

			var contentText = '<div>';
			angular.forEach(content, function (component) {

				// todo
				// does it make more sense to have an ng-repeat for each element

			});
			contentText += '</div>';



			//check for controller, must be already included (e.g. by mixin)
			if (page.controller) {
				var locals = {};
				locals.$scope = $scope;
				var ctrl = $controller(page.controller, locals);
				contentText.data('ngControllerController', ctrl);
			}

			$scope.content = $compile(contentText)($scope);

			(!!page.onload) && console.log('loading onload script');
			return Application.script(page.onload);

		}, function(error) {
			//content loading error
			console.log(error);
			//todo - do correctly
			return $q.reject(load.error(error));
		});
	};





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



	var favorite = function(id) {
		console.log("favorite trail with id: " + id);
	};

	// PERSISTENCE

	var persist = {};

	return {
		loadPage : loadPage,
		compile : compile,
		share : Clotho.share,
		favorite : favorite,
		persist : persist
	}
}]);

Application.Trails.directive('TrailPageComponent', ['$compile', function($compile) {

	var pageComponentTypes = {};

	pageComponentTypes.hint = function loadHint(hint) {
		if (!hint) return $q.when();

		var hintDiv = '<div class="pull-right" hint-button="'+hint+'"></div>';

		return $q.when(hintDiv);
	};

	pageComponentTypes.text = function loadText(text) {
		if (!text) return $q.when();
		return $q.when('<div>' + text + '</div>');
	};

	pageComponentTypes.markdown = function loadMarkdown (text) {
		if (!text) return $q.when();
		return $q.when('<ui-markdown>' + text + '</ui-markdown>');
	};

	pageComponentTypes.video = function loadVideo(obj) {
		if (!obj) return $q.when();

		var videoId = YoutubeService.extract_youtube( (angular.isString(obj) ? obj : obj.id) );
		$scope.videoParams = (!!obj.params) ? obj.params : {};

		console.log(obj, obj.autoplay === true, ((obj.autoplay === true) ? 'false' : (obj.mini || true)));

		//note - need single outer parent div to compile properly after replace (maybe not in ng-1.2.x)
		//add this attr to move to next automatically on-complete="next()"
		var template = '<div><div youtube="' + videoId + '" params="videoParams" start-mini="'+ ((obj.autoplay === true) ? false : (obj.mini || true)) +'" autoplay="'+ (obj.autoplay || true) +'"></div></div>';

		//todo - write to avoid timeout? video doesn't update on next() otherwise - probably need to defer instantiation till later
		return $timeout(function() {
			return template;
		});
	};

	pageComponentTypes.template = function loadTemplate (url) {
		if (!url) return $q.when();

		return $http.get(url, {cache:$templateCache}).then(function (success) {
			return success.data;
		}, function(error) {
			console.log('error retrieving template at: ' + url);
			return '<p class="alert alert-error">That template couldn\'t be found :(</p>'
		});
	};

	pageComponentTypes.quiz = function loadQuiz (content) {
		if (!content) return $q.when();

		// todo - copy this, so don't delete dictionary on original object below

		$scope.quiz = content;

		//extend the scope with Page.quiz.dictionary (separate from Page.dictionary)
		if (content.dictionary) {
			var dict = angular.copy(content.dictionary);
			//deletes reference too, i.e. scope.quiz.dictioary
			delete content.dictionary;
			angular.extend($scope.quiz, dict);
		}

		$scope.gradeCallback = function(data) {
			console.log('quiz grade callback result: ' + data);
		};

		var template = '<div trail-quiz ng-model="quiz" grade-callback="gradeCallback()" advance="next()"></div>';
		return $q.when(template);
	};

	pageComponentTypes.error = function loadError(error) {
		return '<h4>Something didn&apos;t work - that type of Page wasn&apos;t recognized</h4>' +
			error ? '<div>'+error+'</div>' : '';
	};

	var loadPageComponent = function TrailLoadPageComponent(obj, type) {
		if (!type || !angular.isString(type)){
			console.log('no type passed');
			return;
		}

		if (!pageComponentTypes[type]) {
			return pageComponentTypes.error();
		}

		return pageComponentTypes[type](obj);
	};

	return {
		restrict: 'EA',
		scope: false, //avoid isolate
		compile: function compile(tElement, tAttrs, transclude) {
			return {
				pre: function preLink(scope, element, attrs) {
					scope.component = attrs.component;
					scope.componentType = attrs.component.type;
					scope.componentParams = attrs.component.params;

					pageComponentTypes[scope.componentType](scope.componentParams).then(function(result) {
						element.html($compile(result)(scope))
					});
				},
				post: function postLink(scope, element, attrs) {

				}
			}
		}
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

		//don't activate already active one
		if ($scope.current == indices) return;

		console.log('1 - activating Page');

		$scope.current = indices;

		//in form <Chapter>-<Page>
		var pos = indices.split("-");
		var page = $scope.trail.contents[pos[0]]['pages'][pos[1]];

		$scope.content = Trails.loadPage(page, $scope);
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

	//todo ? move this heavy logic into Trails Service

	$scope.next = function() {
		var oldpos = (typeof $scope.current != 'undefined') ? $scope.current.split("-") : [0, -1],
			newpos;

		if (typeof $scope.trail.contents[oldpos[0]]['pages'][+oldpos[1] + 1] != 'undefined')
			newpos = oldpos[0] + '-' + (+oldpos[1] + 1);
		else if (typeof $scope.trail.contents[+oldpos[0] + 1]['pages'] != 'undefined')
			newpos = (+oldpos[0] + 1) + '-' + 0;
		else {
			$scope.current = undefined;
			return;
		}

		$scope.activate(newpos);
	};

	$scope.prev = function() {
		console.log($scope.current);
		if ($scope.current == '0-0') return;

		var oldpos = (typeof $scope.current != 'undefined') ? $scope.current.split("-") : [0, 1],
			newpos;

		if (typeof $scope.trail.contents[oldpos[0]]['pages'][+oldpos[1] - 1] != 'undefined')
			newpos = oldpos[0] + '-' + (+oldpos[1] - 1);
		else if (typeof $scope.trail.contents[+oldpos[0] - 1]['pages'])
			newpos = (+oldpos[0] - 1) + '-' + ($scope.trail.contents[+oldpos[0] - 1]['pages'].length - 1);
		else {
			$scope.current = undefined;
			return;
		}

		$scope.activate(newpos);
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
			autoplay: '@?',
			startMini: '@?',
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
						videoId : scope.videoId,
						playerVars : {
							autoplay : (!scope.autoplay && !scope.startMini) ? 0 : 1,
							autohide : 1,
							rel : 0
						},
						events : {}
					};
					scope.params = angular.extend(defaults, scope.params);

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
//todo - need trail in proper format to test

angular.module('clotho.trails').directive('trailPageComponent', function ($compile, $q, $timeout, $http, $templateCache, YoutubeService) {

	return {
		restrict: 'EA',
		scope: false, //avoid isolate
		compile: function compile(tElement, tAttrs, transclude) {
			return {
				pre: function preLink(scope, element, attrs) {
					//DELEGATION TO TYPE

					var pageComponentTypes = {};

					pageComponentTypes.hint = function loadHint(hint) {
						if (!hint || !angular.isObject(hint)) return $q.when();

						var hintDiv = '<div class="pull-right" hint-button="' + hint + '"></div>';

						return $q.when(hintDiv);
					};

					//note - you can pass in angular bindings since they will be compiled later...
					pageComponentTypes.text = function loadText(text) {
						if (!text) return $q.when();
						return $q.when('<div>' + text + '</div>');
					};

					pageComponentTypes.markdown = function loadMarkdown(text) {
						if (!text) return $q.when();
						return $q.when('<ui-markdown>' + text + '</ui-markdown>');
					};

					pageComponentTypes.video = function loadVideo(obj) {
						if (!obj) return $q.when();

						var videoId = YoutubeService.extract_youtube((angular.isString(obj) ? obj : obj.id));
						scope.videoParams = (!!obj.params) ? obj.params : {};

						scope.videoParams.autoplay = angular.isDefined(obj.autoplay) ? obj.autoplay : false;
						scope.videoParams.mini = angular.isDefined(obj.mini) ? obj.mini : false;

						//note - outer div so compiles properly
						var template = '<div><div youtube="' + videoId + '" params="videoParams"></div></div>';

						//note - may need to introduce $timeout if want to call next() on onComplete from controller
						return $q.when(template);
					};

					pageComponentTypes.template = function loadTemplate(url) {
						if (!url || !angular.isString(url)) return $q.when();
						console.log(url);
						return $q.when('<div ng-include="\'' + url + '\'"></div>');

						/*return $http.get(url, {cache:$templateCache}).then(function (success) {
						 return success.data;
						 }, function(error) {
						 console.log('error retrieving template at: ' + url);
						 return '<p class="alert alert-error">That template couldn\'t be found :(</p>'
						 });*/
					};

					pageComponentTypes.quiz = function loadQuiz(content) {
						if (!content || !angular.isObject(content)) return $q.when();

						scope.quiz = angular.copy(content);

						//extend the scope with Page.quiz.dictionary (separate from Page.dictionary)
						if (scope.quiz.dictionary) {
							var dict = angular.copy(scope.quiz.dictionary);
							delete scope.quiz.dictionary;
							angular.extend(scope.quiz, dict);
						}

						//testing, should also run next() clause on complete
						scope.gradeCallback = function (data) {
							console.log('quiz grade callback result: ' + data);
						};

						var template = '<div trail-quiz ng-model="quiz" grade-callback="gradeCallback()" advance="next()"></div>';
						return $q.when(template);
					};

					pageComponentTypes.error = function loadError(error) {
						return $q.when('<h4>Something didn&apos;t work - that type of Page wasn&apos;t recognized</h4>' +
							error ? '<div>' + error + '</div>' : '');
					};

					var loadPageComponent = function TrailLoadPageComponent(obj, type) {
						if (!type || !angular.isString(type)) {
							console.log('no type passed');
							return;
						}

						if (!pageComponentTypes[type]) {
							return pageComponentTypes.error();
						}

						return pageComponentTypes[type](obj);
					};

					scope.createPageComponent = function() {
						//console.log('creating component (' +  scope.componentType + ')', scope.componentParams);
						return loadPageComponent(scope.componentParams, scope.componentType).then(function (result) {
							return $q.when(element.html($compile(result)(scope)));
						});
					};

					//LOAD IT
					//note- not linked in pre-phase
					//return scope.createPage();
				},
				post: function postLink(scope, element, attrs) {
					scope.$watch(attrs.trailPageComponent, function (oldComp, newComp) {
						if (!!newComp) {
							scope.component = newComp;
							scope.componentType = scope.component.type;
							scope.componentParams = scope.component.params;
							scope.createPageComponent();
						}
					});
				}
			}
		}
	}
});
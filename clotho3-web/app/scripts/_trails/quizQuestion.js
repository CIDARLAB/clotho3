angular.module('clotho.quiz')
/**
 * @ngdoc directive
 * @name quiz-question
 *
 * @description
 * Given a properly formed Quiz Question, displays it using the proper template (chosen by quiz.type)
 *
 * @attr ngModel {Quiz Question}
 * @attr gradeCallback {Function} on attr as: `grade-callback="myCallback($result)"`
 */
	.directive('quizQuestion', function (QuizQuestion, $interpolate, $q) {

		return {
			restrict: "E",
			require: 'ngModel',
			templateUrl : 'views/_trails/quiz/_container.html',
			scope : {
				quiz : '=ngModel',
				gradeCallback : '&?'
			},
			link: function quizQuestionLink(scope, element, attrs) {

				/* setup + utilities */

				var defaultOptions = {
					showAnswer : false,
					allowMultiple : false,
					allowRetry : true,
					randomization : true
				};

				function createEmptyAnswer() {
					var type = scope.quiz.question.type;
					switch (type) {
						case 'mc' : {
							return new Array(scope.quiz.question.options.length);
						}
						default : {
							return  ''
						}
					}
				}

				var interpolateArguments = function (args) {
					return QuizQuestion.interpolateArguments(args, scope);
				};

				function regenerateDynamic () {
					return QuizQuestion.interpolateDictionary(scope.quiz.dictionary)
					.then(function (interpolatedDict) {
						//extend scope with the dictionary
						angular.extend(scope, interpolatedDict);

						//note - options are interpolated on submission, not here (they are displayed as interpolated using directive below
					});
				}

				scope.inputEmpty = function () {
					return angular.isUndefined(scope.$meta) || angular.isEmpty(scope.$meta.input);
				};

				/* quiz grading + retrying etc. */

				scope.grade = function () {
					//can only interpolate strings -- don't want to pass numbers or booleans through
					var interpolatedInput = angular.isString(scope.$meta.input) ? $interpolate(scope.$meta.input)(scope) : scope.$meta.input;
					var interpolatedArgs = interpolateArguments(scope.quiz.grade.args);

					$q.all({
						result : QuizQuestion.grade(scope.quiz, interpolatedInput, interpolatedArgs),
						feedback : QuizQuestion.feedback(scope.quiz, interpolatedInput)
					})
					.then(function (results) {

						//don't want to overwrite the input
						angular.extend(scope.$meta , {
							submitted : true,
							response : !!results.result,
							currentFeedback : results.feedback
						});

						if (angular.isDefined(attrs.gradeCallback)) {
							scope.gradeCallback({
								$input : interpolatedInput,
								$feedback : results.feedback,
								$result: results.result
							});
						}
					});
				};

				scope.reset = function () {
					scope.$meta = {
						input : '',
						submitted : false,
						response : null,
						currentFeedback : null,
						loading : false
					};
				};

				scope.retry = function () {
					regenerateDynamic().then(function () {
						scope.reset();
					});
				};

				//future - make this snazzier...
				scope.showAnswer = function () {
					var interpolatedArgs = interpolateArguments(scope.quiz.grade.args);

					QuizQuestion.grade(scope.quiz, null, interpolatedArgs, true)
					.then(function (answer) {
						scope.$meta.input = answer;
					});
				};

				/* watchers */

				scope.$watch('quiz', function (newquiz) {

					if (!angular.isEmpty(newquiz)) {
						//add default options, in new variable so don't alter quiz itself
						scope.quizOptions = angular.extend({}, defaultOptions, newquiz.options || {});
						scope.retry();
					} else {
						scope.reset();
						scope.$meta.loading = true;
					}
				});
			}
		}
	})
/*
 * internal directives used in quiz questions
 */

//note - this just changes what is shown. Must be interpolated on submission
//functions as ng-bind kinda... extra intpolation later
.directive('qqValue', function ($compile) {
	return {
		link: function (scope, element, attrs) {
			scope.$watch(attrs.qqValue, function (value) {
				var el = angular.element('<span>' + value + '</span>');
				element.html($compile(el)(scope));
			});
		}
	}
})
.directive('qqQuestion', function ($compile) {
	return {
		restrict: "E",
		link: function postLink(scope, element, attrs) {
			//cannot simply compile because the text itself needs to be compiled, so watch for question to change
			scope.$watch('quiz.question.question', function (newval) {
				var wrapped = angular.element('<div>' + (newval || '') + '</div>');
				element.html($compile(wrapped)(scope));
			});
		}
	}
})
.directive('qqTemplate', function ($http, $compile) {
		return {
		restrict: "E",
		link: function (scope, element, attrs) {
			//because we are getting a template, need to set up a watch
			//can't use ng-include becuase sets up isolate scope
			scope.$watch('quiz.question.type', function (newval) {
				if (newval) {
					$http.get('views/_trails/quiz/' + scope.quiz.question.type + '.html', {cache : true})
					.success(function (data, headers) {
						element.html($compile(data)(scope));
					})
					.error(function (data, headers) {
						element.html('<p>template not found</p>');
					});
				}
			});
		}
	}
})
.directive('qqHint', function () {
	return {
		restrict: "E",
		replace: true,
		template: '<div popover="{{quiz.question.hint}}"' +
			'popover-trigger="mouseenter"' +
			'popover-placement="left">' +
			'<span class="glyphicon glyphicon-info-sign"></span>' +
			'</div>'
	}
})
.directive('qqActions', function () {
	return {
		restrict: "E",
		replace: true,
		templateUrl: 'views/_trails/quiz/_actions.html'
	}
});
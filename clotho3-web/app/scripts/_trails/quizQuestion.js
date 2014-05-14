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
	.directive('quizQuestion', function (QuizQuestion, $interpolate) {

		return {
			restrict: "E",
			require: 'ngModel',
			templateUrl : 'views/_trails/quiz/_container.html',
			scope : {
				quiz : '=ngModel',
				gradeCallback : '&?'
			},
			compile: function compile(tElement, tAttrs, transclude) {
				return {
					pre: function preLink(scope, element, attrs) {

						//get template based on type

					},
					post: function postLink(scope, element, attrs) {

						/* setup + utilities */

						//todo - default options

						//todo - hide if not interpolated

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

						//given array / object / *, return array for arguments for function
						function interpolateArguments (args) {
							var interpolated = [];
							angular.forEach(args, function (arg) {
								interpolated.push($interpolate(arg)(scope))
							});
							return interpolated;
						}

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

							QuizQuestion.grade(scope.quiz, interpolatedInput, interpolatedArgs).then(function (result) {

								console.log('result is ' , result);

								var relevantFeedback = QuizQuestion.feedback(scope.quiz, interpolatedInput);

								//don't want to overwrite the input
								angular.extend(scope.$meta , {
									submitted : true,
									response : !!result,
									currentFeedback : relevantFeedback
								});

								if (angular.isDefined(attrs.gradeCallback)) {
									scope.gradeCallback({
										$input : interpolatedInput,
										$feedback : relevantFeedback,
										$result: result
									});
								}
							});
						};

						scope.reset = function () {
							scope.$meta = {
								input : '',
								submitted : false,
								response : null,
								currentFeedback : null
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

						//todo

						//init pending watches
						scope.retry();
					}
				}
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
			//simply compile and binding will update accordingly
			var wrapped = angular.element('<div>' + scope.quiz.question.question + '</div>');
			element.html($compile(wrapped)(scope));
		}
	}
})
.directive('qqTemplate', function ($http, $compile) {
		return {
		restrict: "E",
		link: function (scope, element, attrs) {
			//because we are getting a template, need to set up a watch
			//can't use ng-include becuase sets up isolate scope
			//scope.$watch('quiz.question.type', function (newval) {
			//	console.log('quiz type ' + newval);
				$http.get('views/_trails/quiz/' + scope.quiz.question.type + '.html', {cache : true})
				.success(function (data, headers) {
					element.html($compile(data)(scope));
				})
				.error(function (data, headers) {
					element.html('<p>template not found</p>');
				});
			//});
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
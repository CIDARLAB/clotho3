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
	.directive('quizQuestion', function (QuizQuestion) {

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

						//pull in dictionary
						//interpolate dictionary.dynamic, then question, options, etc.

						/* todo - implement options:

							 options : {
								 checkAnswer : false,
								 showAnswer : false,
								 allowMultiple : false,
								 allowRetry : true,
								 randomization : false
							 }
						 */

						/* setup + utilities */

						scope.createEmptyAnswer = function (quiz, value) {
							//todo - for quizzes which have array as answer, check quiz type
						};

						//todo - create scope for answers etc. outside quiz

						scope.answerUndefined = angular.isEmpty;


						function interpolateArguments (args) {
							if (!angular.isArray(args)) {
								return [];
							}

						}

						/* quiz grading + retrying etc. */

						scope.grade = function () {
							//need to interpolate arguments
							var interpolatedArgs = interpolateArguments(scope.grade.args);
							QuizQuestion.grade(scope.quiz.id, scope.input).then(function (result) {

								//todo

								scope.response = !!result;

								if (angular.isDefined(attrs.gradeCallback)) {
									scope.gradeCallback({$result: result});
								}
							});
						};

						scope.feedback = function () {
							QuizQuestion.feedback()
						};

						scope.reset = function () {

						};

						scope.retry = function () {
							QuizQuestion.retry().then(function (result) {
								//todo
							});
						};

						scope.showAnswer = function () {

						};

						scope.checkAnswer = function () {

						};

						/* watchers */

						//todo - watch for things to interpolate - can watch quiz directly


						//init pending watches
						QuizQuestion.interpolateDictionary(scope.quiz.dictionary)
						.then(function (interpolated) {
							angular.extend(scope, interpolated);
						});


					}
				}
			}
		}
	})

	.directive('quizQuestionQuestion', function ($compile) {
		return {
			restrict: "E",
			scope: false, //do not create isolate scope
			link: function (scope, element, attrs) {
				var wrapped = angular.element('<div>' + scope.quiz.question.question + '</div>');
				element.html($compile(wrapped)(scope));
			}
		}
	})
/**
 * Internal directive to handle the input of the question, for which a specific template is pulled for each question type
 */
.directive('quizQuestionInput', function ($http) {
		return {
			restrict: "E",
			scope: false, //do not create isolate scope
			link: function (scope, element, attrs) {
				//todo - get proper template, tie to answer
			}
		}
	});
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
			scope : {
				quiz : '=ngModel',
				gradeCallback : '&?'
			},
			compile: function compile(tElement, tAttrs, transclude) {
				return {
					pre: function preLink(scope, element, attrs) {

						//pull in dictionary

						//get template based on type

						//interpolate everything (question, answers, etc.)

					},
					post: function postLink(scope, element, attrs) {

						/* setup + utilities */

						scope.createEmptyAnswer = function (quiz, value) {
							//todo - for quizzes which have array as answer, check quiz type
						};

						scope.answerUndefined = angular.isEmpty;

						/* quiz grading + retrying etc. */

						scope.grade = function (quiz) {
							QuizQuestion.grade().then(function (result) {

								//todo

								if (angular.isDefined(attrs.gradeCallback)) {
									scope.gradeCallback({$result: result});
								}
							});
						};

						scope.reset = function () {

						};

						scope.retry = function () {
							QuizQuestion.retry().then(function (result) {
								//todo
							});
						};

						/* watchers */

						//todo - watch for things to interpolate

					}
				}
			}
		}
	});
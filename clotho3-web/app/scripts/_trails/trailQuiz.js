angular.module('clotho.quiz')
/**
 * @ngdoc directive
 * @name trail-quiz
 *
 * @attr ngModel
 * @attr grade-callback {Function}
 * @attr advance {Function}
 */
.directive('trailQuiz', function($http, $templateCache, $compile, Clotho, $interpolate, $q, $sce) {

	//temporary, to get out of API
	//todo - reference by ID
	/**
	 * @name gradeQuiz
	 *
	 * @param {*} questionValue
	 * @param {*} input
	 * @param {string} answerGen ID of function to run to generate answer
	 *
	 * @description
	 * wrapper for grade quiz funciton on server (easier to change in one place)
	 *
	 */
	var gradeQuiz = function (questionValue, input, answerGen) {
		return Clotho.run('gradeQuiz', [questionValue, input, answerGen]);
	};

	return {
		restrict: "EA",
		require: 'ngModel',
		scope: {
			quiz: '=ngModel',
			gradeCallback : '&?',
			advance : '&?' //todo - events?
		},

		compile: function compile(tElement, tAttrs, transclude) {
			return {
				pre: function preLink(scope, element, attrs) {

					//todo -- rethink extending scope with whole quiz (i.e. want dictionary, but maintain quiz namespace?)
					angular.extend(scope, scope.quiz);

					// note - see also grade and retry functions, and load.quiz

					//$watch in postLink
					scope.interpolatedQuestion = $interpolate(scope.quiz.question)(scope);

					scope.parsedQuestion = function() {
						return $sce.trustAsHtml('<h5>' + scope.quiz.question + '</h5>');
					};

					$http.get('views/_trails/quiz/' + scope.quiz.type + '-partial.html', {cache: $templateCache})
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
						gradeQuiz(quiz.questionValue, quiz.answer, quiz.answerGenerator).then(function (result) {
							console.log('gradeQuiz result: ' + result);
							scope.quiz.submitted = true;
							scope.quiz.response = {};
							scope.quiz.response.result = result;

							if (angular.isDefined(attrs.gradeCallback)) {
								scope.gradeCallback({$result : result});
							}
						});
					};

					scope.$watch('quiz.questionValue', function () {
						scope.interpolatedQuestion = $interpolate(scope.quiz.question)(scope)
					});

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
});
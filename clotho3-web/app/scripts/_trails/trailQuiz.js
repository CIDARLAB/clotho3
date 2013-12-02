angular.module('clotho.trails').directive('trailQuiz', function($http, $templateCache, $compile, Clotho, $interpolate, $q) {

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
					// todo - use $parse quiz question
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
});
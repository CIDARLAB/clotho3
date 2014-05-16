angular.module('clotho.quiz')
/**
 * @ngdoc service
 * @name QuizQuestion
 *
 * @description
 * Service to handle Quiz Question grading, retrying, feedback etc.
 *
 * todo - after writing, move to server
 */
	.service('QuizQuestion', function (Clotho, $q, $interpolate, $rootScope) {

		/**
		 * @description
		 * given array or object, interpolate each using passed scope ($scope or object)
		 *
		 * @param args {Array|Object} with values to interpolate. Will not be interpolated if not a string.
		 * @param scope {Object|Scope} dictionary for interpolating
		 * @returns {Array} interpolated array / object (same as input)
		 */
		function interpolateArguments (args, scope) {
			var interpolated = angular.isArray(args) ? [] : {};
			angular.forEach(args, function (arg, key) {
				interpolated[key] = angular.isString(arg) ? $interpolate(arg)(scope) : arg;
			});
			return interpolated;
		}

		this.interpolateArguments = interpolateArguments;

		/**
		 * @description
		 an author can specify either:

		 an answer (grade.answer) - a string, array, or static function
		 a function (grade.function) handling all grading, returning a boolean

		 * @param quiz {QuizQuestion}
		 * @param input {*}
		 * @param args {Array|Object}
		 * @param returnAnswer {boolean}
		 */
		this.grade = function gradeQuestion(quiz, input, args, returnAnswer) {
			// (on server) synchronous get of quiz
			// var quiz = Clotho.get(questionId);

			var result;

			args = _.isArray(args) ? args : [];

			if (quiz.grade.answer) {

				var answerValue = quiz.grade.answer.value;

				if (quiz.grade.answer.type == 'number') {
					if (returnAnswer) {
						return $q.when(answerValue);
					}
					var tolerance = quiz.grade.answer.tolerance;

					result = (input > (answerValue - tolerance) && input < (answerValue + tolerance) );
				}
				else if (quiz.grade.answer.type == 'boolean') {
					if (returnAnswer) {
						return $q.when(answerValue);
					}

					result = (input == !!answerValue);
				}
				if (quiz.grade.answer.type == 'string') {
					if (returnAnswer) {
						return $q.when(answerValue);
					}

					result = input.toLowerCase() == answerValue.toLowerCase();
				}
				else if (quiz.grade.answer.type == 'array') {
					if (returnAnswer) {
						return $q.when(answerValue[0]);
					}

					var foundIndex = _.find(answerValue, function (val) {
						return val.toLowerCase() == input.toLowerCase();
					});
					result = foundIndex >= 0;
				}
				else if (quiz.grade.answer.type == 'function') {
					if (returnAnswer) {
						return Clotho.run(answerValue, args);
					}

					result = Clotho.run(answerValue, args).then(function (fnResult) {
						if (_.isString(fnResult)) {
							return input.toLowerCase() == fnResult.toLowerCase();
						}
						else {
							return input == fnResult;
						}
					});
				}
				else {
					//shouldn't happen in this definition, but may add features later
					Clotho.say('answer was provided in wrong format (cannot grade)');
				}
			}
			else if (quiz.grade.function) {
				if (returnAnswer) {
					return Clotho.run(quiz.grade.function, [input, args]);
				}

				result = Clotho.run(quiz.grade.function, [input, args]);
			}
			else {
				// catch for answer and function undefined (should never happen)
				Clotho.say('no answer was defined...');
				result = null;
			}

			// after determining whether answer is correct
			// future - do additional processing with result
			// persistResponse(questionUUID, input, result);
			// communicateWithClient(questionUUID, result, input, args);

			//finally, return the result
			//return result;

			return $q.when(result);
		};

		//future - handle dynamic feedback
		this.feedback = function generateFeedback (quiz, input) {
			//todo - handle case-insensitive
			if (!quiz.feedback) {
				return $q.when(null);
			}
			if (_.isUndefined(input)) {
				return $q.when(quiz.feedback.default || null);
			}
			var staticFeedback = quiz.feedback.static || {};
			return $q.when(staticFeedback[input] || quiz.feedback.default);
		};


		this.interpolateDictionary = function interpolateDictionary (dictionary) {

			if (angular.isEmpty(dictionary)) {
				return $q.when();
			}

			var interpolatedDict = {};

			//first add static values
			if (angular.isDefined(dictionary.static)) {
				angular.extend(interpolatedDict, dictionary.static);
			}

			//for storing sequential runs
			var prevPromise = $q.when(interpolatedDict);

			//go through dynamic values, sequentially, so previous interpolations can be used
			if ( angular.isDefined(dictionary.dynamic) ) {
				angular.forEach(dictionary.dynamic, function (keyObj) {
					//no good way of just getting the value of a single object
					//don't want to alter the dictionary itself, this creates a copy
					_.mapValues(keyObj, function (runObj, key) {

						// interpolate, submit, add to chain to resolve sequentially
						prevPromise = prevPromise.then(function (intermediateDict) {
							var interpolatedArgs = interpolateArguments(runObj.args, intermediateDict);

							return Clotho.run(runObj.id, interpolatedArgs)
							.then(function (result) {
								//extend interpolated dictionary with new value
								intermediateDict[key] = result;
								return intermediateDict;
							});
						});
					});
				});
			}

			return prevPromise;
		};

		//note - retrying is handled directly in the quiz question since it is only a matter of reinterpolating the dictionary, and extending the scope
	});
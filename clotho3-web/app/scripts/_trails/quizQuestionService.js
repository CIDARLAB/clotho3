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
	.service('QuizQuestion', function (Clotho, $q, $interpolate) {

		/**
		 * @description
		 an author can specify either:

		 an answer (grade.answer) - a string, array, or static function
		 a function (grade.function) handling all grading, returning a boolean

		 * @param questionId
		 * @param input
		 * @param args {Array|Object}
		 */
		this.grade = function gradeQuestion(questionId, input, args) {
			// synchronous get of quiz
			//var quiz = Clotho.get(questionId);

			var result;

			Clotho.get(questionId, {mute: true})
			.then(function (quiz) {

				args = _.isArray(args) ? args : [];

				//todo - make case insensitive

				// check if grade.answer is defined
				if (quiz.grade.answer) {
					// Author supplied a single answer
					if (quiz.grade.answer.type == 'string') {
						//check equality, lowercasing
						result = (input == quiz.grade.answer.value);
					}
					// Author supplied a few possible answers
					else if (quiz.grade.answer.type == 'array') {
						//check index, lowercasing
						result = _.indexOf(quiz.grade.answer.value, input) >= 0;
					}
					// Author supplied a function to run
					else if (quiz.grade.answer.type == 'function') {
						//check against function result
						var functionResult = Clotho.run(quiz.grade.answer.value, args);
						result = (input == functionResult);
					}
					else {
						//shouldn't happen in this definition, but may add features later
					}
				}
				//use grade.function is answer not defined
				else if (quiz.grade.function) {
					result = Clotho.run(quiz.grade.function, args);
				}
				else {
					// catch for answer and function undefined (should never happen)
					Clotho.say('no answer was defined...');
					result = null;
				}
			});

			// after determining whether answer is correct
			// future - do additional processing with result
			// persistResponse(questionUUID, input, result);
			// communicateWithClient(questionUUID, result, input, args);

			//finally, return the result
			//return result;

			return $q.when(result).then(function (det) {
				return det;
			});
		};


		this.retry = function (questionId) {
			//todo
		};


		this.feedback = function (questionId, input, args) {
			//todo
		};


		this.interpolateDictionary = function interpolateDictionary (dictionary) {

			//for storing sequential runs
			var promiseChain = $q.when({});

			if (angular.isEmpty(dictionary)) {
				return promiseChain;
			}

			var interpolatedDict = {};

			//first add static values
			if (angular.isDefined(dictionary.static)) {
				angular.forEach(dictionary.static, function (obj) {
					angular.extend(interpolatedDict, obj);
				});
			}

			//go through dynamic values
			if ( angular.isDefined(dictionary.dynamic) ) {
				angular.forEach(dictionary.dynamic, function (obj) {
					//no good way of just getting the value of a single object
					//don't want to alter the dictionary itself, this creates a copy
					angular.extend(interpolatedDict, _.mapValues(obj, function (val) {
						var interpolatedValue = $interpolate(val)(interpolatedDict);
						console.log(interpolatedValue);

						//todo - check if just needs to interpolate, or also submit
						//once interpolate, need to submit, add to chain of promises to resolve sequentially

						return interpolatedValue;
					}));
				});
			}

			return promiseChain.then(function () {
				return interpolatedDict;
			});
		}
	})
;
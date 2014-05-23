angular.module('clotho.construction')
	.service('ConstructionSimulator', function (ConstructionReactions, $q) {

		/**
		 * @name parseName
		 * @description
		 * Given a file and name, checks for name in dictionary, or retrieves first result from Clotho if not present
		 *
		 * @param file
		 * @param name
		 * @param dictionary
		 *
		 * @returns {Promise} promise, resolving to value
		 */
		function parseName (file, name, dictionary) {
			//first, try the dictionary if its given, return result if find it
			if (!angular.isEmpty(dictionary)) {

			}

			//if we get here, check clotho by query for name
		}

		/**
		 * @name processStep
		 * @description
		 * Compute the value for a given step, looking up the input and running the reaction
		 *
		 * @param file
		 * @param step
		 *
		 * @returns {Promise} promise, resolving to object { <outputName> : <value> }
		 */
		function processStep (file, step) {

		}

		var process = function processConstructionFile (file) {

			//initial checks

			//make sure dictionary is an object
			file.dictionary = angular.isDefined(file.dictionary) ? file.dictionary : {};

			//iterate through steps
			_.forEach(file.steps, function (step) {

				//parse names and add to dictionary
				_.forEach(step.input, function (input) {

				})

			});

		};

		return {
			reactions : reactions,
			process : process
		}

	});
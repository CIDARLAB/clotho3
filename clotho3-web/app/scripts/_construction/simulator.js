angular.module('clotho.construction')
	.service('ConstructionSimulator', function (ConstructionReactions, $q, Debug, Clotho) {

		var Debugger = new Debug('ConstructionSimulator', '#88bb88');

		/**
		 * @name parseName
		 * @description
		 * Given a file and name, checks for name in dictionary, or retrieves first result from Clotho if not present
		 *
		 * @param file {ConstructionFile}
		 * @param value {string}
		 *
		 * @returns {Promise} promise, resolving to value
		 */
		function parseName (file, value) {

			//make sure it's a string
			if (!_.isString(value)) {
				return $q.when(value);
			}

			//try the dictionary if its given, return result if find it
			if (!_.isEmpty(file.dictionary)) {
				var dictionaryResult = file.dictionary[value];
				if (dictionaryResult) {
					return $q.when(dictionaryResult);
				}
			}

			//if we get here, check clotho by query for name, return first result
			return Clotho.query({name : value}, {mute : true}).then(function (results) {
				if (results.length != 1) {
					Debugger.warn('ambiguous query! ' + results.length + ' results for ' + value);
				}
				return results[0];
			});
		}

		/**
		 * @name processStep
		 * @description
		 * Compute the value for a given step, looking up the input and running the reaction
		 *
		 * @param file {ConstructionFile}
		 * @param step {string}
		 *
		 * @returns {Promise} promise, resolving to object { <outputName> : <value> }
		 */
		function processStep (file, step) {

		}

		/**
		 * @name process
		 * @description
		 * Given a well-formed construction file, iterate through each step and calculate final product, adding intermediates to dictionary. returns a copy of the construction file with extended dictionary, and result under the key 'final'
		 *
		 * @param inputFile {ConstructionFile}
		 *
		 * @returns {Promise} promise, resolving to object { <outputName> : <value> }
		 */
		var process = function processConstructionFile (inputFile) {

			/* config */
			var final_key = 'final';

			//create copy, with defaults
			var file = _.assign({
				dictionary : {}
			}, inputFile);

			//initial checks, return if nothing to do
			if (_.isEmpty(file.steps)) {
				return file;
			}

			/* simulate */

			//iterate through steps
			_.forEach(file.steps, function (step) {

				//parse names and add to dictionary
				_.forEach(step.input, function (input) {

					//todo - handle nested (i.e. array)

				});

			});

			//create final key
			var lastKey = _.last(file.steps).output;
			if (lastKey != final_key) {
				file.dictionary[final_key] = file.dictionary[lastKey];
			}

			return file;

		};

		return {
			process : process
		}

	});
angular.module('clotho.construction')
	.service('ConstructionSimulator', function (ConstructionReactions, $q, Debug, Clotho) {

		var Debugger = new Debug('ConstructionSimulator', '#ee8888');

		/* config */
		var FINAL_DICT_KEY = 'final';

		function disambiguateOnClotho (value) {
			//make sure it's a string - just return booleans and numbers
			if (!_.isString(value)) {
				return $q.when(value);
			}

			//check clotho by query for name, return first result
			return Clotho.query({name : value}, {mute : true}).then(function (results) {
				if (results.length != 1) {
					Debugger.log('ambiguous query! ' + results.length + ' results for ' + value);
				}
				return results[0];
			});
		}

		/**
		 * @name parseName
		 * @description
		 * Given a file and name, checks for name in dictionary, or retrieves first result from Clotho if not present
		 *
		 * Does NOT handle arrays
		 *
		 * @param file {ConstructionFile}
		 * @param value {string}
		 *
		 * @returns {Promise} promise, resolving to value
		 */
		function parseName (file, value) {

			//make sure it's a string - just return booleans and numbers - don't want to look up in dictionary
			if (!_.isString(value)) {
				return $q.when(value);
			}

			//try the dictionary if its given, return result if find it
			if (!_.isEmpty(file.dictionary)) {
				var dictionaryResult = file.dictionary[value];
				if (!_.isUndefined(dictionaryResult)) {
					return $q.when(dictionaryResult);
				}
			}

			//if we get here, check clotho by query for name, return first result
			return disambiguateOnClotho(value);
		}

		/**
		 * @name processStepArgs
		 *
		 * @description
		 * Given a step, flattens input and and returns map of input name to promised value
		 *
		 * @param file {ConstructionFile}
		 * @param step
		 * @returns {object} Object where keys are flatted input, values are values (promises)
		 */
		function processStepArgs (file, step) {

			//determine all arguments
			var inputArgs = _.flatten(step.input);

			//retrieve from dictionary / database, saving as array
			var promiseArray = _.map(inputArgs, function (arg) {
				return parseName(file, arg);
			});

			return _.zipObject(inputArgs, promiseArray)
		}

		/**
		 * Given a step and populated dictionary, interpolate the names of step input
		 * note - for now, we're just going one level in
		 */
		function constructStepArgs (file, step) {

			//helper function to handle primitive / string / array of args
			//arg should is either a primitive (just return it), string (look it up), or array (recurse)
			function interpolateArg (array) {
				return _.map(array, function (arg) {
					return _.isString(arg) ? file.dictionary[arg] :
								 _.isArray(arg) ? interpolateArg(arg) : arg;
				});
			}

			return interpolateArg(step.input);
		}

		/**
		 * @name processStep
		 * @description
		 * Compute the value for a given step, looking up the input and running the reaction
		 *
		 * Assumes that dictionary is already populated for all arguments needed in this step (by running processStepArgs)
		 *
		 * @param file {ConstructionFile}
		 * @param step
		 *
		 * @returns {Promise} promise, resolving to object { <outputName> : <value> }
		 */
		function processStep (file, step) {

			var reaction = ConstructionReactions[step.reaction];

			//everything already exists in the dictionary, so we just need to pull those values
			var parsedInput = constructStepArgs(file, step);

			return Clotho.run(reaction.reactionId, parsedInput, {mute : true})
			.then(function(result) {
				var obj = {};
				obj[step.output] = result;
				return obj;
			});
		}

		/**
		 * @name preprocess
		 * @description
		 * Given a well-formed construction file, iterate through each step and interpolate arguments.
		 *
		 * @param inputFile {ConstructionFile}
		 * @param createNew {boolean} If true, passed file will be copied and not edited itself
		 *
		 * @returns {Promise} promise, resolving to object { <outputName> : <value> }
		 *
		 */
		//todo - check all inputs, deferring to dictionary
		var interpolateInputs = function (inputFile, createNew) {

			//create copy, with defaults
			var file = _.assign({
				dictionary : {}
			}, !!createNew ? _.clone(inputFile, true) : inputFile);

			//initial checks, return if nothing to do
			if (_.isEmpty(file.steps)) {
				return $q.when(file);
			}

			var dictionaryTerms = _.keys(file.dictionary);
			var clothoInputs = {};

			Debugger.groupCollapsed('Processing Input Arguments');

			_.forEach(file.steps, function (step) {
				dictionaryTerms.push(step.output);

				_.forEach(_.flatten(step.input), function (input) {
					if ( _.isString(input) && _.indexOf(dictionaryTerms, input) < 0 ) {
						clothoInputs[input] = disambiguateOnClotho(input);
					}
				});
			});

			return $q.all(clothoInputs)
			.then(function (resolvedPartialDict) {
				Debugger.groupEnd();
				_.assign(file.dictionary, resolvedPartialDict);
				return file;
			});
		};

		/**
		 * @name process
		 * @description
		 * Given a well-formed construction file, iterate through each step and calculate final product, adding intermediates to dictionary.
		 *
		 * @param inputFile {ConstructionFile}
		 * @param createNew {boolean} If true, passed file will be copied and not edited itself
		 *
		 * @returns {Promise} promise, resolving to object { <outputName> : <value> }
		 *
		 */
		var process = function processConstructionFile (inputFile, createNew) {

			//create copy, with defaults
			var file = _.assign({
				dictionary : {}
			}, !!createNew ? _.clone(inputFile, true) : inputFile);

			//initial checks, return if nothing to do
			if (_.isEmpty(file.steps)) {
				return $q.when(file);
			}

			//note - on client this is async, so need promises
			var prevStepPromise = $q.when();

			/* simulate */

			//group debugger
			Debugger.groupCollapsed('Simulating File');

			//iterate through steps, sequentially - need to handle sequential promises (like getting schema promises)
			_.forEach(file.steps, function (step, index) {

				//add this step to the promise chain
				prevStepPromise = prevStepPromise.then(function() {

					Debugger.log('step ' + index, file.dictionary, step.input, step.output);

					//resolve all input arguments
					return $q.all(processStepArgs(file, step))
				})
				.then(function (resolved) {

					Debugger.log('step ' + index + ' args resolved', resolved);

					//add to dictionary if non-primitive
					_.forEach(resolved, function (arg, key) {
						if ( (_.isArray(arg) || _.isObject(arg)    &&
									_.isUndefined(file.dictionary[key]))  )  {
							file.dictionary[key] = arg;
						}
					});

					//run the step
					return processStep(file, step);
				})
				.then(function (stepResult) {
					Debugger.log('step ' + index + ' processed', stepResult);
					//add the step result to the dictionary
					return $q.when(_.assign(file.dictionary, stepResult));
				});
			});

			return prevStepPromise
			.then(function (chain) {
				//once we've hit this, the promise chain has resolved and file is computed

				//ungroup debugger
				Debugger.groupEnd();

				//create final key
				var lastKey = _.last(file.steps).output;
				if (lastKey != FINAL_DICT_KEY) {
					file.dictionary[FINAL_DICT_KEY] = file.dictionary[lastKey]
				}

				//hack - remove cut marks
				if (_.isArray(file.dictionary[FINAL_DICT_KEY])) {
					file.dictionary[FINAL_DICT_KEY][0].sequence = file.dictionary[FINAL_DICT_KEY][0].sequence.replace(/[\^\|_]/ig, '');
				} else if (_.isObject(file.dictionary[FINAL_DICT_KEY]) && !_.isUndefined(file.dictionary[FINAL_DICT_KEY].sequence)) {
					file.dictionary[FINAL_DICT_KEY].sequence = file.dictionary[FINAL_DICT_KEY].sequence.replace(/[\^\|_]/ig, '');
				}

				return file;
			});
		};

		return {
			interpolateInputs : interpolateInputs,
			process : process
		}

	});
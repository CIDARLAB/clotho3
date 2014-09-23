angular.module('clotho.commandbar')
.service('CommandBar', function(Clotho, ClientAPI, ClothoCommandHistory, Debug, ClothoSchemas, $timeout, $q) {

	/******* config ******/
	var options = {
		dateFilter : 'mediumTime',
		timeFilter : 'timestamp'
	};

	var Debugger = new Debug('Command Bar', "#ffbb55");

	/******* elements ******/


	var getCommandBarElement = function() {
		return angular.element($document[0].querySelector('[clotho-command-bar]'));
	};

	var getTokenizerElement = function () {
		return angular.element($document[0].querySelector('[clotho-command-bar] [clotho-tokenizer]'));
	};

	//todo - should capture from commandBar directive as possible, e.g. via controller
	//note - call as needed, ensure exists in DOM
	var getCommandBarInput = function () {
		return angular.element($document[0].querySelector('[clotho-command-bar] [clotho-reference-autocomplete]'));
	};

	var commandBarInputModel = 'query';

	var getCommandBarLogButton = function () {
		return angular.element($document[0].getElementById('clotho_logButton'));
	};

	var focusInput = function() {
		getCommandBarInput().focus();
	};

	var setInput = function (string) {
		getCommandBarInput().scope().setQueryString(string);
	};


	/****** display ******/
	var display = {};
	display.query = '';
	display.log = false; // activity log
	display.logSnippet = false; // snippet right of log button

	display.toggle = function(field, value) {
		display[field] = angular.isDefined(value) ? value : !display[field];
	};

	display.toggleActivityLog = function () {
		display.log = !display.log;
		if (display.log) {
			log.unread = '';
		}
	};

	//submitting is high up so children can access

	var submit = function (input) {
		if (angular.isEmpty(input) || !angular.isObject(input)) {
			input = display.query;
		}

		//remove trailing whitespace
		input.query = angular.isDefined(input.query) ? input.query.trim() : '';

		if (!!input.query) {
			var submission = {
				class : 'info',
				from : 'client',
				text: input.query,
				timestamp : Date.now()
			};

			/*
			//DEPRECATED - will allow element to check for delimiter as described in GH#399
			//check if we have a sharable
			if (angular.isDefined(input.tokens)) {
				submission.tokens = angular.map(input.tokens, function (token) {
					return token.model;
				});
			}
			*/

			ClientAPI.say(submission);

			return Clotho.submit(input).then(function(result){
				display.query = '';
				ClientAPI.say({text: result, class: 'success'});
			}, function (rejection) {
				//do nothing...
			});
		}
		else {
			return $q.when(false);
		}
	};

	/****** facade ******/

	return {
		options : options,
		display : display,
		setQuery : function(item, $event) {
			if (angular.isDefined( $event )) {
				$event.preventDefault();
			}

			if (!item) return;

			display.query = !angular.isEmpty(item.value) ? item.value : item.text;
		},
		submit : submit,

		//interaction
		getCommandBarElement: getCommandBarElement,
		getTokenizerElement : getTokenizerElement,
		getCommandBarInput : getCommandBarInput,
		commandBarInputModel : commandBarInputModel,
		focusInput : focusInput,
		setInput : setInput
	}

});

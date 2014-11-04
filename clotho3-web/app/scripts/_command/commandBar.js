angular.module('clotho.commandbar')
.service('CommandBar', function(Clotho, ClientAPI, ClothoCommandHistory, Debug, ClothoSchemas, $timeout, $q, $document) {

	var Debugger = new Debug('Command Bar', "#ffbb55");

	/******* elements ******/


	var getCommandBarElement = function() {
		return angular.element($document[0].querySelector('clotho-command-bar'));
	};

	var getTokenizerElement = function () {
		return angular.element($document[0].querySelector('clotho-command-bar [clotho-tokenizer]'));
	};

	//todo - should capture from commandBar directive as possible, e.g. via controller
	var getCommandBarInput = function () {
		return angular.element($document[0].querySelector('clotho-command-bar [clotho-reference-autocomplete]'));
	};

	var commandBarInputModel = 'query';

	var getCommandBarLogButton = function () {
		return angular.element($document[0].getElementById('clotho_logButton'));
	};

	var focusInput = function() {
		getCommandBarInput().focus();
	};

  //brittle...
	var setInput = function (string) {
		getCommandBarInput().scope().setQueryString(string);
	};

	/****** display ******/
	var display = {};
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

	/****** facade ******/

	return {
		display : display,
		setQuery : setInput,

		//interaction
		getCommandBarElement: getCommandBarElement,
		getTokenizerElement : getTokenizerElement,
		getCommandBarInput : getCommandBarInput,
		commandBarInputModel : commandBarInputModel,
		focusInput : focusInput,
		setInput : setInput
	}

});

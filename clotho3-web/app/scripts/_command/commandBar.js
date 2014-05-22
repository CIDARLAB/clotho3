angular.module('clotho.commandbar')
.service('CommandBar', function(Clotho, ClientAPI, Debug, $timeout, $q, $document) {

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

	//todo - should capture from commandBar directive as possible, e.g. via controller
	var getCommandBarInput = function () {
		return angular.element($document[0].querySelector('[clotho-command-bar] [clotho-autocomplete]'));
	};

	var getCommandBarLogButton = function () {
		return angular.element($document[0].getElementById('clotho_logButton'));
	};

	var focusInput = function() {
		getCommandBarInput().focus();
	};


	/******* log data *******/
	var log = {};

	log.entries = [
		{
			"text" : "Welcome to Clotho!",
			"from" : "server",
			"class" : "success",
			"timestamp" : Date.now()
		}
	];

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

	/***** functions *****/

	log.timeout = null;

	function receiveMessage (data) {
		log.unread = (!!log.unread && !display.log) ? log.unread + 1 : 1;
		log.entries.unshift(data);
		Debugger.log('LOG - entries: ', log.entries);
		display.toggle('logSnippet', true);
		log.startLogTimeout();
	}

	log.startLogTimeout = function() {
		log.cancelLogTimeout();

		log.timeout = $timeout( function() {
			display.toggle('logSnippet', false);
		}, 10000);
	};

	log.cancelLogTimeout = function() {
		$timeout.cancel(log.timeout);
	};

	var submit = function (input) {
		console.log('input is "' + input + '"');
		if (angular.isEmpty(input) || !angular.isObject(input)) {
			input = display.query || '';
		}

		//remove trailing whitespace
		input.query = input.query.trim();

		/*
		 Debugger.log(query);
		 Debugger.log(display.queryHistory);
		 Debugger.log(log.entries)
		 */

		if (!!input.query) {
			var submission = {class : 'info', from : 'client', text: input.query, timestamp : Date.now()};

			ClientAPI.say(submission);

			return Clotho.submit(input).then(function(result){
				console.log('resetting query');
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

	/****** listeners *****/

	//called by clientAPI.say(), so will be called indirectly in submit a few lines above
	Clotho.listen("activityLog", function (data) {
		receiveMessage(data);
	}, 'searchbar');


	/****** facade ******/

	return {
		options : options,
		display : display,
		log : log,
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
		getCommandBarInput : getCommandBarInput,
		focusInput : focusInput
	}

});
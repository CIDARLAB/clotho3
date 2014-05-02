//rename directive

angular.module('clotho.commandbar').service('CommandBar', function(Clotho, ClientAPI, Debug, $timeout, $q, $document) {

	/******* config ******/
	var options = {
		dateFilter : 'mediumTime',
		timeFilter : 'timestamp'
	};

	var Debugger = new Debug('Command Bar', "#ffbb55");

	/******* elements ******/

	//todo - should capture from commandBar directive as possible

	var getCommandBarElement = function() {
		return angular.element($document[0].querySelector('[clotho-command-bar]'));
	};

	var getCommandBarInput = function () {
		return angular.element($document[0].getElementById('clotho_command_input'));
	};

	var getCommandBarLogButton = function () {
		return angular.element($document[0].getElementById('clotho_logButton'));
	};

	var focusInput = function() {
		getCommandBarInput().focus();
	};

	function showActivityLog () {
		display.log = true;
	}



	/******* log data *******/
	var log = {};

	var autocomplete = {};
	autocomplete.autocompletions = [];
	autocomplete.autoDetail = {};
	autocomplete.detailTemplate = {};
	autocomplete.detailModel = {};
	autocomplete.detailUUID = -1;

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
	display.queryHistory = [];
	display.autocomplete = false; // autocomplete list
	display.autocompleteDetail = false; //pane to left of autocomplete
	display.autocompleteDetailInfo = false; // e.g. command or author
	display.help = false; // help menu far right
	display.log = false; // activity log
	display.logSnippet = false; // snippet right of log button

	//basic toggles
	display.show = function (field) {
		if (!display[field])
			display[field] = true;
	};

	display.hide = function(field) {
		if (display[field])
			display[field] = false;
	};

	display.toggle = function(field) {
		display[field] = !display[field];
	};

	// todo - should be CSS
	display.genLogPos = function() {
		var target = getCommandBarLogButton()[0];
		display.logpos = {
			left : (target.offsetLeft + (target.scrollWidth / 2) - 200) + "px",
			top : (target.offsetTop + target.scrollHeight)  + "px"
		};
	};

	/*
	//note - now in css
	display.genAutocompletePos = function() {
		var target = getCommandBarInput()[0];
		display.autocompletePos = {
			left : (target.offsetLeft) + "px",
			top : (target.offsetTop + target.clientHeight)  + "px"
		};
	};
	*/


	display.detail = function(uuid) {
		if (typeof uuid == 'undefined') return;

		if (uuid != autocomplete.detailUUID) {
			display.hide('autocompleteDetailInfo');
			autocomplete.detailUUID = uuid;
		}

		Clotho.autocompleteDetail(autocomplete.detailUUID).then(function(result) {
			autocomplete.autoDetail = result;
			display.show('autocompleteDetail');
		});
	};

	display.undetail = function() {
		display.hide('autocompleteDetail');
		display.hide('autocompleteDetailInfo');
		autocomplete.detailModel = {};
	};

	//todo - avoid using index in case sort - have to namespace
	display.detailInfo = function (type, index) {
		//choose template
		switch (type) {
			case 'command' : {
				autocomplete.detailTemplate = 'views/_command/detail-command.html';
				break;
			}
			case 'author' : {
				autocomplete.detailTemplate = 'views/_command/detail-author.html';
				break;
			}
			default : {}
		}
		//choose model
		autocomplete.detailModel = autocomplete.autoDetail.sharables[index];
		if (type == "author")
			autocomplete.detailModel = autocomplete.detailModel.author;

		display.show('autocompleteDetailInfo');
	};


	/***** functions *****/

	log.timeout = null;

	function receiveMessage (data) {
		log.unread = (!!log.unread && !display.log) ? log.unread + 1 : 1;
		log.entries.unshift(data);
		Debugger.log('LOG - entries: ', log.entries);
		display.show('logSnippet');
		log.startLogTimeout();
	}

	log.startLogTimeout = function() {
		log.cancelLogTimeout();

		log.timeout = $timeout( function() {
			display.hide('logSnippet');
		}, 10000);
	};

	log.cancelLogTimeout = function() {
		$timeout.cancel(log.timeout);
	};

	var execute = function (command) {
		Debugger.log("execute called (not implemented) " + command);
		display.hide('autocomplete');
		display.undetail();
	};

	var submit = function (query) {
		if (!query) {
			query = display.query;
		}

		/*
		 Debugger.log(query);
		 Debugger.log(display.queryHistory);
		 Debugger.log(log.entries)
		 */

		if (!!query) {
			var submission = {class : 'info', from : 'client', text: query, timestamp : Date.now()};
			display.queryHistory.push(submission);

			//display.autocomplete = false;
			display.undetail();

			ClientAPI.say(submission);

			//note - temporary, patch as object pending tokenizer
			query = {
				query : query,
				tokens : []
			};

			return Clotho.submit(query).then(function(result){
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
			if (typeof $event != 'undefined')
				$event.preventDefault();
			/*if (item.type != 'command') {
			 display.undetail();
			 }*/
			if (!item) return;

			display.query = !!item.value ? item.value : item.text;
		},
		autocomplete : autocomplete,
		submit : submit,
		execute : execute,

		//interaction
		getCommandBarElement: getCommandBarElement,
		getCommandBarInput : getCommandBarInput,
		focusInput : focusInput,
		showActivityLog : showActivityLog
	}

});
angular.module('clotho.interface').service('$focus', function($document, $timeout, $q, $parse, Clotho, CommandBar) {

	//note - this is a function
	var searchBarInput = CommandBar.getCommandBarInput;

	var maxZ = function(selector) {
		return Math.max(0, Math.max.apply(null, angular.map($document[0].querySelectorAll(selector || "*"),
			function (v) {
				return parseFloat(angular.element(v).css("z-index")) || null;
			})
		));
		/*
		// jquery version (not in use)
		return Math.max.apply(null,
			$.map($('body *'), function(e,n) {
				if (angular.element(e).css('position') != 'static')
					return parseInt(angular.element(e).css('z-index')) || 1;
			})
		);
		*/
	};

	var setZ = function(zindex, element) {
		return $q.when(element.css({"z-index": zindex, 'position' : 'relative'}));
	};

	var bringToFront = function(element) {
		var oldZ = element.css('z-index');
		var newZ = maxZ() + 1;
		setZ(newZ, element);
		return $q.when(oldZ);
	};

  //note - if pass model, assumes it is on element's $scope
	var typeOut = function(element, string, model) {
		var inputsVal = {input: true, textarea : true, select: true},
			timeOut,
			scope,
			valType = (!!inputsVal[angular.lowercase(element[0].nodeName)]) ? "val" : "text",
			txtLen = string.length,
			charInd = 0,
			curString,
			deferred = $q.defer();

		if (!!model) {
			scope = element.scope();
		}

		function typeIt() {
			timeOut = $timeout(function() {
				charInd++;
				curString = string.substring(0, charInd);

				if (scope) {
					scope.$apply($parse(model).assign(scope, curString));
				} else {
					element[valType](curString + '|');
				}

				typeIt();

				if (charInd == txtLen) {
					if (scope) {
						//we're fine
					} else {
						element[valType](element[valType]().slice(0, -1)); // remove the '|'
					}

					deferred.resolve();
					$timeout.cancel(timeOut);
				}

			}, Math.round(Math.random() * 50) + 30);
		}

		typeIt();
		return deferred.promise;
	};

	var typeOutSearch = function(string, submit) {

		//create single element for this function
		var commandInput = searchBarInput();

		string = angular.isArray(string) ? string : [string];

    commandInput.focus();

    return typeOut(commandInput, string, CommandBar.commandBarInputModel)
      .then(function () {
        //clean up reference
        commandInput.remove();
      });
	};

	var backdrop = angular.element("<div>").addClass('modal-backdrop fade');
	//gross
	backdrop.bind('click', function (e) {
		e.preventDefault();
		removeBackdrop();
	});

	//this is gross and hacky
	var addBackdrop = function(zindex) {
		$document.find('body').append(backdrop.css("z-index", zindex || maxZ() + 1));
		return $timeout(function() {backdrop.addClass('in')});
	};

	var removeBackdrop = function() {
		return $q.when(backdrop.removeClass('in'))
			.then(function() {
				return $timeout(function() {
					backdrop.remove()
				}, 150)
			});
	};

	//return function to un-highlight
	var highlightElement = function(el) {
		var oldZ = el.css('z-index');

		addBackdrop();
		setZ(maxZ() + 1, el);

		return $q.when(function() {
			setZ(oldZ, el);
			removeBackdrop();
		});
	};

	//future- move to angular animation
	var pulseElement = function(el) {
		var deferred = $q.defer();
		el.animate( { backgroundColor: "ffff99" }, 300 )
			.animate( { backgroundColor: "transparent" }, 300, function() {
				deferred.resolve();
			});

		return deferred.promise;
	};

	return {
		maxZ : maxZ,
		setZ : setZ,
		bringToFront : bringToFront,
		typeOut : typeOut,
		typeOutSearch : typeOutSearch,
		addBackdrop : addBackdrop,
		removeBackdrop : removeBackdrop,
		highlightElement : highlightElement,
    pulseElement : pulseElement
	}
});

//note - jQuery reliance
//todo - rewrite for autonomy
angular.module('clotho.interface').service('$focus', function($document, $timeout, $q, Clotho) {

	var searchBarInput = ($('#clothoCommandBarInput'));

	var maxZ = function() {
		return Math.max.apply(null,
			$.map($('body *'), function(e,n) {
				if ($(e).css('position') != 'static')
					return parseInt($(e).css('z-index')) || 1;
			})
		);
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



	var typeOut = function(element, string, model) {
		var inputsVal = {input: true, textarea : true, select: true},
			valType = (!!inputsVal[angular.lowercase(element[0].nodeName)]) ? "val" : "text",
			timeOut,
			txtLen = string.length,
			charInd = 0,
			deferred = $q.defer();


		function setDescendentProperty(obj, desc, val) {
			var arr = desc.split(".");
			while(arr.length > 1 && (obj = obj[arr.shift()]));
			obj[arr.shift()] = val;
		}

		function typeIt() {
			timeOut = $timeout(function() {
				charInd++;
				element[valType](string.substring(0, charInd) + '|');
				typeIt();

				if (charInd == txtLen) {
					element[valType](element[valType]().slice(0, -1)); // remove the '|'

					//update scope
					if (!!model) {
						var scope = element.scope();
						//todo - handle two layers in
						setDescendentProperty(scope, model, string);
						scope.$apply();
					}

					deferred.resolve();
					$timeout.cancel(timeOut);
				}

			}, Math.round(Math.random() * (30 - 30)) + 30);
		}

		typeIt();
		return deferred.promise;
	};

	var typeOutSearch = function(string, submit) {

		return $q.when(searchBarInput.focus())
			.then(function() {
				return highlightElement(searchBarInput)
			})
			.then(function(unhighlight) {

				//todo - handle array of strings to input

				return typeOut(searchBarInput, string, 'display.query')
					.then(function() {
						return unhighlight;
					});
			})
			.then(function(unhighlight) {
				return $timeout(function() {unhighlight()}, 600).then(function() {
					if (submit) {
						//CommandBar.display.log = true;
						//return CommandBar.submit(string);
						$q.when(searchBarInput.parents('form').submit())
					} else {
						return $q.when(searchBarInput.focus())
					}
				});
			});
	};


	var backdrop = angular.element("<div>").addClass('modal-backdrop fade');

	//this is gross and hacky
	var addBackdrop = function(zindex) {

		/*backdrop.bind('click', function (e) {
		 e.preventDefault();
		 removeBackdrop();
		 });*/

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
		el
			.animate( { backgroundColor: "ffff99" }, 300 )
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
		highlightElement : highlightElement
	}
});
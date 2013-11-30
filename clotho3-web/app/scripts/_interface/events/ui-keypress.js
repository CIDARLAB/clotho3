/**
 @NOTE - not distributed via package manager, need to update manually:
 https://github.com/angular-ui/ui-utils/blob/master/modules/keypress/keypress.js
 @lastUpdated 11/18/13
 @note custom keypress service at end

 * @description Bind one or more handlers to particular keys or their combination
 * @param hash {mixed} keyBindings Can be an object or string where keybinding expression of keys or keys combinations and AngularJS Exspressions are set.
 *
 * Object syntax: "{ keys1: expression1 [, keys2: expression2 [ , ... ]]}".
 * String syntax: ""expression1 on keys1 [ and expression2 on keys2 [ and ... ]]"".
 *
 * Expression is an AngularJS Expression, and key(s) are dash-separated combinations of keys and modifiers (one or many, if any. Order does not matter).
 * Supported modifiers are 'ctrl', 'shift', 'alt' and key can be used either via its keyCode (13 for Return) or name.
 * Named keys are 'backspace', 'tab', 'enter', 'esc', 'space', 'pageup', 'pagedown', 'end', 'home', 'left', 'up', 'right', 'down', 'insert', 'delete', 'period', 'comma'.
 *
 * @example

 <input ui-keypress="{enter:'x = 1', 'ctrl-shift-space':'foo()', 'shift-13':'bar()'}" />
 <input ui-keypress="foo = 2 on ctrl-13 and bar('hello') on shift-esc" />

 * @note keypress for arrows and some other keys will not work -- use keydown or keyup in that case
 **/

angular.module('ui.keypress',[]).
	factory('keypressHelper', function keypress($parse, $document){
		var keysByCode = {
			8: 'backspace',
			9: 'tab',
			13: 'enter',
			27: 'esc',
			32: 'space',
			33: 'pageup',
			34: 'pagedown',
			35: 'end',
			36: 'home',
			37: 'left',
			38: 'up',
			39: 'right',
			40: 'down',
			45: 'insert',
			46: 'delete'
		};

		var capitaliseFirstLetter = function (string) {
			return string.charAt(0).toUpperCase() + string.slice(1);
		};

		return function(mode, scope, elm, attrs) {
			var params, combinations = [];
			params = scope.$eval(attrs['ui'+capitaliseFirstLetter(mode)]);

			//CUSTOM
			if (elm == $document) {
				params = attrs;
			}

			// Prepare combinations for simple checking
			angular.forEach(params, function (v, k) {
				var combination, expression;
				expression = $parse(v);

				angular.forEach(k.split(' '), function(variation) {
					combination = {
						expression: expression,
						keys: {}
					};
					angular.forEach(variation.split('-'), function (value) {
						combination.keys[value] = true;
					});
					combinations.push(combination);
				});
			});

			// Check only matching of pressed keys one of the conditions
			//custom - assign to handler so can also use in $scope.$destroy
			var handler = function (event) {
				// No need to do that inside the cycle
				var metaPressed = !!(event.metaKey && !event.ctrlKey);
				var altPressed = !!event.altKey;
				var ctrlPressed = !!event.ctrlKey;
				var shiftPressed = !!event.shiftKey;
				var keyCode = event.keyCode;

				// normalize keycodes
				if (mode === 'keypress' && !shiftPressed && keyCode >= 97 && keyCode <= 122) {
					keyCode = keyCode - 32;
				}

				// Iterate over prepared combinations
				angular.forEach(combinations, function (combination) {

					var mainKeyPressed = combination.keys[keysByCode[keyCode]] || combination.keys[keyCode.toString()];

					var metaRequired = !!combination.keys.meta;
					var altRequired = !!combination.keys.alt;
					var ctrlRequired = !!combination.keys.ctrl;
					var shiftRequired = !!combination.keys.shift;

					if (
						mainKeyPressed &&
							( metaRequired === metaPressed ) &&
							( altRequired === altPressed ) &&
							( ctrlRequired === ctrlPressed ) &&
							( shiftRequired === shiftPressed )
						) {
						// Run the function
						scope.$apply(function () {
							combination.expression(scope, { '$event': event });
						});
					}
				});
			};

			elm.bind(mode, handler);

			//custom
			//kill on scope desctruction, if pass in scope
			scope.$on('$destroy', function() {
				elm.unbind(mode, handler);
			});
		};
	});

/**
 * Bind one or more handlers to particular keys or their combination
 * @param hash {mixed} keyBindings Can be an object or string where keybinding expression of keys or keys combinations and AngularJS Exspressions are set. Object syntax: "{ keys1: expression1 [, keys2: expression2 [ , ... ]]}". String syntax: ""expression1 on keys1 [ and expression2 on keys2 [ and ... ]]"". Expression is an AngularJS Expression, and key(s) are dash-separated combinations of keys and modifiers (one or many, if any. Order does not matter). Supported modifiers are 'ctrl', 'shift', 'alt' and key can be used either via its keyCode (13 for Return) or name. Named keys are 'backspace', 'tab', 'enter', 'esc', 'space', 'pageup', 'pagedown', 'end', 'home', 'left', 'up', 'right', 'down', 'insert', 'delete'.
 * @example <input ui-keypress="{enter:'x = 1', 'ctrl-shift-space':'foo()', 'shift-13':'bar()'}" /> <input ui-keypress="foo = 2 on ctrl-13 and bar('hello') on shift-esc" />
 **/
angular.module('ui.keypress').directive('uiKeydown', function(keypressHelper){
	return {
		link: function (scope, elm, attrs) {
			keypressHelper('keydown', scope, elm, attrs);
		}
	};
});

angular.module('ui.keypress').directive('uiKeypress', function(keypressHelper){
	return {
		link: function (scope, elm, attrs) {
			keypressHelper('keypress', scope, elm, attrs);
		}
	};
});

angular.module('ui.keypress').directive('uiKeyup', function(keypressHelper){
	return {
		link: function (scope, elm, attrs) {
			keypressHelper('keyup', scope, elm, attrs);
		}
	};
});


/**
 * @note CUSTOM
 * @example
 * var x = $keypress.on('keypress', {'enter' : 'foo()'}, $scope);
 * ...
 * $keypress.off(x);
 * @note
 * you can't do keypress with alt, control, option, up, down, left, right -- use keydown or keyup
 */
angular.module('ui.keypress').service('$keypress', function(keypressHelper, $document){
	return {
		on : function(mode, actions, scope) {
			return keypressHelper(mode, scope, $document, actions);
		},
		//note - should be handled automatically on $scope.$destroy (so use controllers), this is for manual use
		//format: [elm, mode, handler, combinations];
		off:function(handle) {
			handle[0].unbind(handle[1], handle[2]);
		}
	}
});
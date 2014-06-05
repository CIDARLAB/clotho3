
angular.module('clotho.core')
/**
 * @name Debug
 *
 * @description
 * debugging service used across app
 *
 * @usage
 * var Debugger = new Debug(<namespace (defaults to anonymous)>, <hex_color (optional)>);
 * ...
 * Debugger.log('heres a log message');
 * Debugger.info('heres some info');
 * Debugger.warn('heres a warning');
 * Debugger.error('heres an error');
 * Debugger.debug('heres some debugging');
 * Debugger.table(myObject, fieldsToPluckAsArray);
 * Debugger.$log //patch through to angular $log service
 * Debugger.console // direct access to console
 *
 * // Will output message under namespace given, with color given
 */
	.factory('Debug', function($log, $window, $filter) {

		//set to true for chrome to call console.trace()
		var tracingEnabled = false;

		//set to false to disable all debugging
		var disableAll = false;
		//add a module name to prevent its writing to the console.
		var disabledNamespaces = [];

		//store all debug messages here
		var messages = {};

		function namepaceEnabled (namespace) {
			if (disableAll === true) {
				return false;
			} else {
				return disabledNamespaces.indexOf(namespace) < 0;
			}
		}

		return function (namespace, color) {

			namespace = angular.isString(namespace) ? namespace : 'ANONYMOUS';
			color = /(^#[0-9A-F]{6}$)|(^#[0-9A-F]{3}$)/i.test(color) ? color : '#'+Math.floor(Math.random()*16777215).toString(16);

			messages[namespace] = [];
			var debugFunctionality = {};

			angular.forEach(['log', 'warn', 'error', 'debug', 'info'], function (term) {
				//first argument should be message, multiple arguments supported and passed through to console
				debugFunctionality[term] = function () {
					
					//add stack trace, slice out this function
					//this will work in chrome, firefox, and node
					var stack = (new Error()).stack.split('\n');
					// Chrome includes a single "Error" line, FF doesn't.
					if (stack[0].indexOf('Error') === 0) {
						stack = stack.slice(1);
					}
					//slice out this wrapper's context
					stack = stack.slice(1);

					var x = new Date();
					var readableDate = $filter('date')(x, 'hh:mm:ss.sss');

					//collect all messages
					messages[namespace].push({
						message : arguments,
						time : x.valueOf()
					});

					//check if emabled, add stack if option passed
					if (namepaceEnabled(namespace)) {
						$log[term].apply(null,
							//add debugger name
							['%c' + readableDate + ' - %O - ' + namespace + '\t',
									'color: '+ color +';',
									stack
							]
							.concat(Array.prototype.slice.call(arguments, 0))
						);
						if (tracingEnabled) {
							console.trace();
						}
					}
				}
			});

			//need to handle 'table' function separately
			angular.forEach(['table'], function (term) {
				debugFunctionality[term] = $window.console[term];
			});

			//easily pretty-print objects
			debugFunctionality.object = function (obj) {
				messages[namespace].push({
					message : obj,
					time : Date.now().valueOf()
				});

				if (angular.isObject(obj)) {
					$log.log('%c' + JSON.stringify(obj, null, 2), 'color: '+ color +';')
				} else {
					$log.log('%c' + obj, 'color: '+ color +';');
				}
			};

			//simple wrapping for grouping input
			debugFunctionality.group = function (term) {
				if (angular.isDefined($window.console.group)) {
					$window.console.group('%c' + namespace + '\t' + term || "Collapsing Output", 'color: '+ color +';');
				}
			};

			debugFunctionality.groupCollapsed = function (term) {
				if (angular.isDefined($window.console.groupCollapsed)) {
					$window.console.groupCollapsed('%c' + namespace + '\t' + term || "Collapsing Output", 'color: ' + color + ';');
				}
			};

			debugFunctionality.groupEnd = function () {
				if (angular.isDefined($window.console.groupEnd)) {
					$window.console.groupEnd();
				}
			};

			//for strings to be interpolated etc.
			debugFunctionality.$log = $log;
			debugFunctionality.console = $window.console || {};

			return debugFunctionality;
		}
	});

angular.module('clotho.core')
/**
 * decorate the angular $log to include the table message for chrome
 */
	.config(function ($provide) {
		$provide.decorator('$log', ['$delegate', function($delegate){
			$delegate.table = function() {
				var args = [].slice.call(arguments);
				if(window.console && window.console.table)
					console.table(args[0], args[1]);
				else
					$delegate.log(null, args)
			};
			return $delegate;
		}]);
	})
/**
 * @name Debug
 *
 * @description
 * debugging service used across app
 *
 * @usage
 * new Debug(<namespace (defaults to anonymous)>, <hex_color (optional)>);
 * ...
 * Debug.log('heres a log message');
 * Debug.info('heres some info');
 * Debug.warn('heres a warning');
 * Debug.error('heres an error');
 * Debug.debug('heres some debugging');
 * Debug.table(myObject, fieldsToPluckAsArray);
 *
 * // Will output message under namespace given, with color given
 */
	.factory('Debug', function($log) {

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
				debugFunctionality[term] = function (msg) {

					var x = new Date();
					var readableDate = x.getHours() + ':' + x.getMinutes() + ':' + x.getSeconds() + '.' + x.getMilliseconds();

					//collect all messages
					messages[namespace].push({
						message : arguments,
						date : readableDate
					});

					//add functionality
					if (namepaceEnabled(namespace)) {
						$log[term].apply(null, ['%c' + readableDate + ' - ' + namespace + '\t', 'color: '+ color +';'].concat(Array.prototype.slice.call(arguments, 0)));
					}
				}
			});

			debugFunctionality.table = function (obj) {
				$log.table(obj);
			};

			debugFunctionality.object = function (obj) {
				if (angular.isObject(obj)) {
					$log.log('%c' + JSON.stringify(obj, null, 2), 'color: '+ color +';')
				} else {
					$log.log('%c' + obj, 'color: '+ color +';');
				}
			};

			return debugFunctionality;
		}
	});
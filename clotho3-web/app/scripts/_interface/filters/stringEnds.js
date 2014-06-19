angular.module('clotho.interface').filter('stringEnds', function() {
	return function (input, length) {

		length = length || 20;
		var process = function(input) {
			return input.substr(0, length) + '...' + input.substring(input.length - length);
		};

		if (!angular.isString (input) ) {
			if (angular.isArray(input)) {
				var first = angular.isString(input[0]) ? process(input[0]) : '';
				return '<Array>' + (first ? ' First Entry:\n' + first : '');
			}
			else if (angular.isObject(input)) {
				var key = angular.isString(input['key']) ? process(input['key']) : '';
				return '<Object>' + (key ? 'Object.key:\n' + key : '');
			}
			else return input
		} else {
			if (input.length <= length * 2) {
				return input;
			}
		}

		return process(input)
	}
});
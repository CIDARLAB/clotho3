angular.module('clotho.interface').filter('joinArray', function() {
	return function (input, joiner) {
		if (!angular.isObject(input)) {
			return input;
		} else if (angular.isArray(input)) {
			return input.join(joiner);
		} else {
			var output = '';
			angular.forEach(input, function (val, key) {
				output += key + joiner;
			});
			return output.substring(0, output.length - joiner.length);
		}
	}
});
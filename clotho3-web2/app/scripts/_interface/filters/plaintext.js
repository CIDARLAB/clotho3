//basic regex to remove HTML tags
angular.module('clotho.interface').filter('plainText', function() {
	// note -- not perfect, but much faster than jQuery.text()
	// http://jsperf.com/regex-replace-vs-jquery-text
	return function(input, filterOn) {
		return (filterOn) ? input.replace(/(<([^>]+)>)/ig, '') : input;
	}
});
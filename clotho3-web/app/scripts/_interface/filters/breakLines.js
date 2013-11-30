//returns a string joined by 'joiner' after 'charCount' characters
angular.module('clotho.interface').filter('breakLines', function() {
	return function(input, charCount, joiner) {
		//todo - ignore HTML
		return (input.match(new RegExp('.{1,'+charCount+'}', 'gi')) || []).join(joiner || '\n');
	}
});
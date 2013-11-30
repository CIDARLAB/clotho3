angular.module('clotho.commandbar').directive('logEntries', function() {

	return {
		restrict: 'A',
		templateUrl: '../views/_command/logEntries.html',
		scope: {
			entries: '=logEntries'
		}
	}
});